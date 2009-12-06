#include <burns.h>
#include <monty.h>
#include <primes.h>
#include <gmp.h>

burns::burns(int bits)
{
	// The number of bits gained per prime is
	// more than 63.999999999999 for the first
	// hundred thousand primes.
	// Take one bit extra to accomodate for this
	// slight deviation from 64.
	uint64 n = (bits + 1) / 64 + 1;

	// Find the primes and construct the Montgomery fields
	fields.reserve(n);
	for(uint64 p = max_uint64; fields.size() < n; p -= 2)
	{
		if(is_prime(p))
		{
			fields.push_back(monty(p));
		}
	}

	// Construct M and M_approx
	mpz_init(M);
	mpz_set_ui(M, 1);
	for(int i=0; i < n; i++)
	{
		uint64 modulus = fields[i].modulus();
		mpz_mul_ui(M, M, modulus);
	}
	// gmp_printf(" M = %Zd\n", M);

	// Construct the (m / M) mod m
	for(int i=0; i < n; i++)
	{
		monty& m = fields[i];
		uint64 M = m.one();
		for(int j = 0; j < n; j++)
		{
			if(i == j) continue;
			uint64 modulus = m.set(fields[j].modulus());
			M = m.mul(M, modulus);
		}
		m.c = m.get(m.get(m.inv(M)));
	}

	Mmod3 = 1;
	for(int i=0; i < n; i++)
	{
		monty& m = fields[i];
		Mmod3 *= m.modulus() % 3;
		Mmod3 %= 3;
	}
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

vector<uint64> burns::set(mpz_t X) const
{
	vector<uint64> xi;
	xi.reserve(size());
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];
		uint64 residue = mpz_fdiv_ui(X, m.modulus());
		xi.push_back(m.set(residue));
	}
	return xi;
}

void burns::get(mpz_t X, vector<uint64> x) const
{
	mpz_t Mi;
	mpz_init(Mi);

	mpz_set_ui(X, 0);
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];

		// Mi = M / m_i
		mpz_divexact_ui(Mi, M, m.modulus());

		uint64 y = m.get_CRT(x[i]);

		// X += Mi y
		mpz_addmul_ui(X, Mi, y);
	}
	mpz_clear(Mi);

	// Reduce mod M
	mpz_mod(X, X, M);
}

uint64 burns::fractional(vector<uint64> x) const
{
	uint64 X = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];
		uint64 Xi = div128(m.get_CRT(x[i]), 0, m.modulus());
		X += Xi;
	}
	return X;
}

uint64 burns::mod64(vector<uint64> x) const
{
	uint64 Xm = 0;
	uint64 Mm = 1;

	uint64 XM = 0;
	uint64 W = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];

		// Xi = X m M⁻¹  mod m
		uint64 Xi = m.get_CRT(x[i]);

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
		const monty& m = fields[i];

		uint64 xi = x[i];

		// xi += (M - (M mod 3)) / 3 mod m
		uint64 Md3 = m.mul(m.sub(m.zero(), m.set(Mmod3)), m.inv(m.set(3)));
		xi = m.add(xi, Md3);

		// Xi = X m M⁻¹  mod m
		uint64 Xi = m.get_CRT(xi);

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

uint64 burns::mod(vector<uint64> x, const monty& k) const
{
	uint64 Xk = k.zero();
	uint64 Mk = k.one();

	uint64 XM = 0;
	uint64 W = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];
		uint64 mk = k.set(m.modulus());

		// Xi = X m M⁻¹  mod m
		uint64 Xi = m.get_CRT(x[i]);

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

	// Try again, but add (M - delta) / 3
	Xk = 0;
	XM = 0;
	W = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];

		uint64 xi = x[i];

		// xi += (M - (M mod 3)) / 3 mod m
		uint64 Md3 = m.mul(m.sub(m.zero(), m.set(Mmod3)), m.inv(m.set(3)));
		xi = m.add(xi, Md3);

		// Xi = X m M⁻¹  mod m
		uint64 Xi = m.get_CRT(xi);

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

/// @returns max_uint64 when the number
uint64 burns::to_uint64(vector<uint64> x) const
{
	uint64 X = fields[0].get(x[0]);
	for(int i=1; i < size(); i++)
	{
		const monty& m = fields[i];
		if(m.get(x[i]) != X)
		{
			return max_uint64;
		}
	}
	return X;
}
