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
	
	// Construct the (m / M) mod m
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
		mrc.push_back(m.get(m.get(m.inv(M))));
	}
	
	// Calculate M mod 3
	Mmod3 = 1;
	for(int i=0; i < n; i++)
	{
		const monty& m = mb.field(i);
		Mmod3 *= m.modulus() % 3;
		Mmod3 %= 3;
	}
	
	// Calculate m_i^-1 mod m_j
	/*
	mji.resize(size());
	for(int j = size() - 1; j >= 0; j--)
	{
		mji[j] = vector<uint64>(j);
		const monty& mj = mb.field(j);
		for(int i=0; i < j; i++)
		{
			const monty& mi = mb.field(i);
			mji[j][i] = mi.inv(mi.get(mj.modulus()));
		}
	}
	*/
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
	
	// TODO: Faster algo using divide and conquer?
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
	// TODO: Faster algo using divide and conquer?
	mpz_t Mi;
	mpz_init(Mi);

	mpz_set_ui(X, 0);
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);

		// Mi = M / m_i
		mpz_divexact_ui(Mi, M, m.modulus());

		uint64 y = m.modular_mul(x[i], mrc[i]);

		// X += Mi y
		mpz_addmul_ui(X, Mi, y);
	}
	mpz_clear(Mi);

	// Reduce mod M
	mpz_mod(X, X, M);
}

void burns::to_mixed_radix(vector<uint64>& x) const
{
	for(int j = size() - 1; j >= 0; j--)
	{
		// const vector<uint64>& mi_inv = mji[j];
		const monty& mj = mb.field(j);
		x[j] = mj.get(x[j]);
		for(int i = 0; i < j; i++)
		{
			
			const monty& mi = mb.field(i);
			x[i] = mi.sub(x[i], mi.get(x[j]));
			// x[i] = mi.mul(x[i], mi_inv[i]);
			
			// The inversion is a slow operation,
			// how can we avoid it without precalculating
			// a O(n^2) sized table?
			
			x[i] = mi.div(x[i], mi.get(mj.modulus()));
		}
	}
	// X = x_1 M/m_1 + x_2 M/m_1m_2 + ... + x_n
}

void burns::from_mixed_radix(vector<uint64>& x) const
{
	for(int i = 0; i < size(); i++)
	{
		const monty& mi = mb.field(i);
		for(int j = 0; j <i; j++)
		{
			const monty& mj = mb.field(j);
			x[j] = mj.mul(x[j], mj.get(mi.modulus()));
			x[j] = mj.add(x[j], mj.get(x[i]));
		}
		x[i] = mi.set(x[i]);
	}
}

int burns::compare(const vector<uint64>& x, const vector<uint64>& y)
{
	// First try fractional
	uint64 xf = fractional(x);
	uint64 yf = fractional(y);
	if(xf > yf) return  1;
	if(xf < yf) return -1;
	
	// Try equality
	if(equals(x, y)) return 0;
	
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

uint64 burns::fractional(const vector<uint64>& x) const
{
	uint64 X = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);
		uint64 Xi = div128(m.modular_mul(x[i], mrc[i]), 0, m.modulus());
		X += Xi;
	}
	return X;
}

uint64 burns::mod64(const vector<uint64>& x) const
{
	uint64 Xm = 0;
	uint64 Mm = 1;
	
	uint64 XM = 0;
	uint64 W = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = mb.field(i);

		// Xi = X m M⁻¹  mod m
		uint64 Xi = m.modular_mul(x[i], mrc[i]);

		uint64 XiM = div128(Xi, 0, m.modulus());
		XM += XiM;
		if(XM < Xi) W++;

		// Mm *= m mod 2⁶⁴
		Mm *= m.modulus();

		// Xm += Xi / m  mod 2⁶⁴
		Xm += Xi * inv_mod64(m.modulus());
	}

	if((XM + size()) > XM)
	{
		// Xm -= W mod 2⁶⁴
		Xm -= W;

		// Xm *= M mod 2⁶⁴
		Xm *= Mm;

		// W is certain, return the result
		return Xm;
	}

	// Try again, but add (M - delta) / 3
	Xm = 0;
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
		uint64 Xi = m.modular_mul(xi, mrc[i]);

		uint64 XiM = div128(Xi, 0, m.modulus());
		XM += XiM;
		if(XM < Xi) W++;

		// Xm += Xi / m  mod 2⁶⁴
		Xm += Xi * inv_mod64(m.modulus());
	}

	// Xm -= W mod 2⁶⁴
	Xm -= W;

	// Xm *= M mod 2⁶⁴
	Xm *= Mm;

	// Xm -= (M - 1) / 3 mod 2⁶⁴
	Xm -= (Mm - Mmod3) * inv_mod64(3);

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
		uint64 Xi = m.modular_mul(x[i], mrc[i]);

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
		uint64 Xi = m.modular_mul(xi, mrc[i]);

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

/// @returns max_uint64 when the number
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
