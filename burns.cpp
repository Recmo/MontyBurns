#include <montyburns.h>
#include <burns.h>
#include <monty.h>
#include <primes.h>
#include <gmp.h>

burns::burns()
{
	nsize = 0;
}

burns::burns(int n)
{
	nsize = n;
	
	// Construct M and M_approx
	mpz_init(M);
	mpz_set_ui(M, 1);
	for(int i=0; i < size(); i++)
	{
		uint64 modulus = mb.field(i).modulus();
		mpz_mul_ui(M, M, modulus);
	}
	// gmp_printf(" M = %Zd\n", M);
	
	// Construct the (m / M) mod m in monty form
	mrc.reserve(size());
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		uint64 M = m.one();
		for(int j = 0; j < n; j++)
		{
			if(i == j) continue;
			uint64 modulus = m.set(mb.field(j).modulus());
			M = m.mul(M, modulus);
		}
		mrc.push_back(m.inv(M));
	}
	
	// Calculate M mod 2^64
	Mmod64 = 1;
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		Mmod64 *= m.modulus();
	}
	
	// Calculate (M + 1) / 2 mod 2^64
	vector<uint64> M2;
	M2.reserve(size());
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		M2[i] = m.inv(m.set_small(2));
	}
	M2mod64 = 0;//mod64(M2); // TODO: Fix it!
}

uint64 div128(uint64 ah, uint64 al, uint64 d)
{
	uint64 q, r;
	asm("divq %4" : "=a"(q), "=d"(r) : "a"(al), "d"(ah), "rm"(d));
	return q;
}

uint64 inv_mod64(uint64 b)
{
	uint64 p = 1;
	for(int i=64; i; i--)
	{
		p *= b;
		b *= b;
	}
	return p;
}

const uint64 ten19 = 10000000000000000000ul;
const uint64 inv3_mod64 = inv_mod64(3);
const uint64 inv10_mod64 = inv_mod64(ten19);
const double inv64 = 5.421010862427522E-20;

vector<uint64> burns::set_uint64(const uint64 X) const
{
	vector<uint64> x;
	x.reserve(size());
	for(int j = 0; j < size(); j++)
	{
		const monty& m = mb.field(j);
		x.push_back(m.set(X));
	}
	return x;
}

vector<uint64> burns::set(const uint64* X, int length) const
{
	if(length == 0)
	{
		return set_uint64(0);
	}
	
	vector<uint64> x;
	x.reserve(size());
	for(int j = 0; j < size(); j++)
	{
		const monty& m = mb.field(j);
		uint64 R2 = m.monty_R() * m.monty_R();
		uint64 xj = m.set(X[length - 1]);
		for(int i = length - 2; i >= 0; --i)
		{
			xj = m.mul(xj, R2);
			xj = m.add(xj, m.set(X[i]));
		}
		x.push_back(xj);
	}
	return x;
}

vector<uint64> burns::set(mpz_t X) const
{
	// int l = X->_mp_size;
	// if (l < 0) l = -l;
	// return set(X->_mp_d, l);
	
	vector<uint64> xi;
	xi.reserve(size());
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		uint64 residue = mpz_fdiv_ui(X, m.modulus());
		xi.push_back(m.set(residue));
	}
	return xi;
}

void burns::get(mpz_t X, const vector<uint64>& x) const
{
	mpz_t Mi;
	mpz_init(Mi);
	
	mpz_set_ui(X, 0);
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		
		// Mi = M / m_i
		mpz_divexact_ui(Mi, M, m.modulus());
		
		uint64 y = m.get(m.mul(x[i], mrc[i]));
		
		// X += Mi y
		mpz_addmul_ui(X, Mi, y);
	}
	mpz_clear(Mi);

	// Reduce mod M
	mpz_mod(X, X, M);
}

/// @brief In place conversion to a mixed radix representation
///
/// Sets x such that X = x_1 M/m_1 + x_2 M/m_1m_2 + ... + x_n
void burns::to_mixed_radix(vector<uint64>& x) const
{
	// TODO: develop a left-to-right method
	// using the algorithm from fractional
	for(int j = size() - 1; j >= 0; j--)
	{
		const monty& m = mb.field(j);
		// Incorporate an extra monty R in factor
		// to eliminate some get and set operations
		uint64 factor = m.neg(m.monty_R() * m.monty_R());
		for(int i = size() - 1; i > j; i--)
		{
			// x += -a_i * m_0 * m_1 * ... m_i-1
			x[j] = m.add(x[j], m.mul(x[i], factor));
			
			// factor *= m_i
			uint64 modi = mb.field(i).modulus();
			factor = m.mul(factor, m.set(modi));
		}
		x[j] = m.div(x[j], m.neg(factor));
	}
}

void burns::from_mixed_radix(vector<uint64>& a) const
{
	for(int i = 0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		uint64 xi = m.set(a[size() - 1]);
		// Incorporate an extra monty R in factor
		// to eliminate some get and set operations
		uint64 factor = m.monty_R() * m.monty_R();
		for(int j = size() - 2; j >= i; j--)
		{
			uint64 modj = mb.field(j + 1).modulus();
			factor = m.mul(factor, m.set(modj));
			xi = m.add(xi, m.mul(a[j], factor));
		}
		a[i] = xi;
	}
}

int burns::compare(const vector<uint64>& x, const uint64 y)
{
	if(fast_false(y >= mb.field(size()-1).modulus()))
	{
		// If the maximum x < y then
		if(size() == 1) return -1;
		
		// Fall back to a full comparisson
		// Mixed radix can be fast since only
		// the two least significants digits are
		// required and a test if the rest is zero
		
		// TODO: Implement
	}
	
	// y < moduli
	uint64 value = mb.field(0).get(x[0]);
	for(int i=1; i < size(); i++)
	{
		const monty& m = mb.field(i);
		if(m.get(x[i]) != value)
		{
			// x >= moduli
			return 1;
		}
	}
	
	// x = value
	if(value > y) return 1;
	if(value < y) return -1;
	return 0;
}

/// Does an unsigned compare
int burns::compare(const vector<uint64>& x, const vector<uint64>& y)
{
	// First try fractional
	uint64 xf = fractional(x);
	uint64 yf = fractional(y);
	if(xf + size() > xf && yf + size() > yf)
	{
		if(xf > yf + size()) return  1;
		if(xf + size() < yf) return -1;
	}
	
	// Try equality
	if(xf == yf && equals(x, y)) return 0;
	
	// Then compare mixed radix
	vector<uint64> xmr = x;
	vector<uint64> ymr = y;
	to_mixed_radix(xmr);
	to_mixed_radix(ymr);
	for(int i = 0; i < size(); i++)
	{
		if(xmr[i] > ymr[i]) return  1;
		if(xmr[i] < ymr[i]) return -1;
	}
	
	// Error
	return 0;
}

uint64 burns::count_wraps(const vector<uint64>& x) const
{
	uint64 Xm = 0;
	uint64 XM = 0;
	uint64 W = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		
		// Xi = X m M⁻¹  mod m
		uint64 Xi = m.get(m.mul(x[i], mrc[i]));
		
		uint64 XiM = div128(Xi, 0, m.modulus());
		XM += XiM;
		if(XM < Xi) W++;
		
		// Xm += Xi / m  mod 2⁶⁴
		Xm += Xi * inv_mod64(m.modulus());
	}
	
	if((XM + size()) > XM)
	{
		// Xm -= W mod 2⁶⁴
		Xm -= W;
		
		// Xm *= M mod 2⁶⁴
		Xm *= Mmod64;
		
		// W is certain, return the result
		return W;
	}
	
	// W is uncertain, could be either W or W + 1
	// Resolve by calculating mod64
	// What are the chances of a mod64 collision?
	
	// Try again, but add (M + 1) / 2
	Xm = 0;
	XM = 0;
	W = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		
		uint64 xi = x[i];
		
		// xi += (M + 1) / 2 mod m = 2^-1 mod m
		// TODO: create a m.half() function, is it faster?
		xi = m.add(xi, m.inv(m.set_small(2)));
		
		// Xi = X m M⁻¹  mod m
		uint64 Xi = m.get(m.mul(xi, mrc[i]));
		
		uint64 XiM = div128(Xi, 0, m.modulus());
		XM += XiM;
		if(XM < Xi) W++;
		
		// Xm += Xi / m  mod 2⁶⁴
		Xm += Xi * inv_mod64(m.modulus());
	}
	
	// Xm -= W mod 2⁶⁴
	Xm -= W;
	
	// Xm *= M mod 2⁶⁴
	Xm *= Mmod64;
	
	// Xm -= (M + 1) / 2 mod 2⁶⁴
	Xm -= M2mod64;
	
	return Xm;
}

/// Returns a number y such that
/// X = (y ... y+size()) * M / 2^64
uint64 burns::fractional(const vector<uint64>& x) const
{
	uint64 X = 0;
	for(int i=0; i < size(); i++)
	{
		// Add (X m M⁻¹  mod m) 2⁶⁴ / m
		// This is not trivially obvious
		// TODO: document
		const monty& m = mb.field(i);
		X += m.mul(x[i], mrc[i]) * m.monty_k();
	}
	return X;
}

uint64 burns::mod64(const vector<uint64>& x) const
{
	// TODO: Fix it!
	
	// Assuming small numbers are more likely
	// then it is faster to start with the
	// shifted version
	uint64 Xm = 0;
	uint64 XM = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		
		uint64 xi = x[i];
		
		// xi += (M + 1) 2⁻¹ mod m = 2⁻¹ mod m
		xi = m.add(xi, m.half());
		
		// XiM = (X m M⁻¹  mod m) 2⁶⁴ / m
		// Xim = (X m M⁻¹  mod m) mod 2⁶⁴
		uint64 Xi = m.mul(xi, mrc[i]);
		uint64 XiM = Xi * m.monty_k();
		uint64 Xim = -m.get(Xi) * m.monty_k();
		
		XM += XiM;
		if(XM < XiM) --Xm;
		Xm += Xim ;
	}
	if((XM + size()) > XM)
	{
		// Xm *= M mod 2⁶⁴
		Xm *= Mmod64;
		
		// Xm -= (M + 1) / 2 mod 2⁶⁴
		Xm -= M2mod64;
		
		// return Xm;
	}
	
	
	cout << endl;
	
	// X is close to (M + 1) / 2 try again
	Xm = 0;
	XM = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		
		// XiM = (X m M⁻¹  mod m) 2⁶⁴ / m
		// Xim = (X m M⁻¹  mod m) mod 2⁶⁴
		uint64 Xi = m.mul(x[i], mrc[i]);
		uint64 XiM = Xi * m.monty_k();
		uint64 Xim = -m.get(Xi) * m.monty_k();
		
		cout << "x[i] = " << m.get(x[i]) << endl;
		cout << "xim  = " << Xim * Mmod64 << endl;
		
		
		XM += XiM;
		if(XM < XiM){ --Xm; cout << "--Xm" << endl; }
		Xm += Xim ;
	}
	
	cout << "XM = " << XM << endl;
	cout << "Xm = " << Xm << endl;
	cout << "Mmod64 = " << Mmod64 << endl;
	
	// Xm *= M mod 2⁶⁴
	Xm *= Mmod64;
	
	// W is certain, return the result
	return Xm;
}

uint64 burns::mod(const vector<uint64>& x, const monty& k) const
{
	uint64 Xk = k.zero();
	uint64 Mk = k.one();
	
	// Try directly calculating the modulus
	uint64 XM = 0;
	uint64 W = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		uint64 mk = k.set(m.modulus());

		// Xi = X m M⁻¹  mod m
		uint64 Xi = m.get(m.mul(x[i], mrc[i]));

		uint64 XiM = div128(Xi, 0, m.modulus());
		XM += XiM;
		if(XM < Xi) W++;

		// Mk *= m mod k
		Mk = k.mul(Mk, mk);

		// Xk += Xi / m  mod k
		Xk = k.add(Xk, k.mul(k.set(Xi), k.inv(mk)));
	}

	if((XM + size()) > XM)
	{
		// Xk -= W mod k
		Xk = k.sub(Xk, k.set(W));

		// Xk *= M mod k
		Xk = k.mul(Xk, Mk);

		// W is certain, return the result
		return Xk;
	}
	
	// The number of wraps was uncertain, x is very close to zero
	//
	// Try again, but add (M - delta) / 3
	Xk = 0;
	XM = 0;
	W = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);

		uint64 xi = x[i];

		// xi += (M - (M mod 3)) / 3 mod m
		uint64 Md3 = m.mul(m.sub(m.zero(), m.set(Mmod3)), m.inv(m.set(3)));
		xi = m.add(xi, Md3);

		// Xi = X m M⁻¹  mod m
		uint64 Xi = m.get(m.mul(xi, mrc[i]));

		uint64 XiM = div128(Xi, 0, m.modulus());
		XM += XiM;
		if(XM < Xi) W++;

		// Xk += Xi / m  mod k
		Xk = k.add(Xk, k.mul(k.set(Xi), k.inv(k.set(m.modulus()))));
	}

	// Xk -= W mod k
	Xk = k.sub(Xk, k.set(W));

	// Xk *= M mod k
	Xk = k.mul(Xk, Mk);

	// Xm -= (M - 1) / 3 mod k
	Xk = k.sub(Xk, k.mul(k.sub(Mk, k.set(Mmod3)), k.inv(k.set(3))));

	return Xk;
}

/// Calculates x = x / y mod M
void burns::safe_div(vector<uint64>& x, const vector<uint64>& y) const
{
	// Ignore montys where y_i = 0
	// Calculate the division
	// Base extend the result to fill in the missing montys
}

/// Calculates x /= y
void burns::safe_div(vector<uint64>& x, const uint64 y, uint64& remainder) const
{
	remainder = mod(x, y);
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		x[i] = m.div(m.sub(x[i], m.set(remainder)), m.set(y));
		// TODO: handle y == m_i
	}
}

/// @returns max_uint64 when X >= 2^64-1
uint64 burns::to_uint64(const vector<uint64>& x) const
{
	uint64 X = mb.field(0).get(x[0]);
	for(int i=1; i < size(); i++)
	{
		const monty& m = mb.field(i);
		if(m.get(x[i]) != X)
		{
			return max_uint64;
		}
	}
	return X;
}

bool burns::equals(const vector<uint64>& x, const uint64 y) const
{
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		if(x[i] != m.set(y))
		{
			return false;
		}
	}
	return true;
}

bool burns::equals_small(const vector<uint64>& x, const uint64 y) const
{
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		if(x[i] != m.set_small(y))
		{
			return false;
		}
	}
	return true;
}

bool burns::equals_zero(const vector<uint64>& x) const
{
	for(int i=0; i < size(); i++)
	{
		if(x[i])
		{
			return false;
		}
	}
	return true;
}

bool burns::equals(const vector<uint64>& x, const vector<uint64>& y) const
{
	for(int i=0; i < size(); i++)
	{
		if(x[i] != y[i])
		{
			return false;
		}
	}
	return true;
}

bool burns::can_shrink(const vector<uint64>& x) const
{
	// y = Shrink x by one monty
	// y = Base extend y by one monty
	// If x == y then return true
	
	// vector<uint64> shrunken = x;
	// shrunken.pop();
	// return shrunken.mod(m[size()]) == x[size()];
	
	
	// Alternative:
	//
	// if(fractional(x) == 0) return true;
}
