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

	// Xm -= (M - 1) / 2 mod 2⁶⁴
	Xm -= (Mm - Mmod3) * inv_mod64(3);

	return Xm;
}

uint64 burns::count_wraps(vector<uint64> x) const
{
	uint64 wraps = 0;
	double X = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];
		double y = m.get_CRT(x[i]);
		X += y / m.modulus();
		if(X >= 1.0)
		{
			wraps++;
			X -= 1.0;
		}
	}
	return wraps;
}

/*
uint64 burns::mod64(vector<uint64> x) const
{
	uint64 wraps = 0;
	uint64 Xm = 0;
	uint64 Mm = 1;
	double X = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& mi = fields[i];

		// Mm *= m_i mod 2⁶⁴
		Mm *= mi.modulus();

		// y_i = x_i × m_i × M^-1  mod m_i
		uint64 y = mi.get_CRT(x[i]);

		// Xm += y / mi  mod 2⁶⁴
		Xm += y * inv64(mi.modulus());

		// X += y / mi  mod 1
		X += static_cast<double>(y) / mi.modulus();
		if(X >= 1.0)
		{
			wraps++;
			X -= 1.0;
		}
	}
	// Xm -= wraps mod m
	Xm -= wraps;

	// Xm *= M mod m
	Xm *= Mm;

	return Xm;
}
*/

uint64 burns::mod(vector<uint64> x, const monty& m) const
{
	uint64 wraps = 0;
	uint64 Xm = m.zero();
	uint64 Mm = m.one();

	double X = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& mi = fields[i];

		// Mm *= m_i mod m
		Mm = m.mul(Mm, m.set(mi.modulus()));

		// y_i = x_i × m_i × M^-1  mod m_i
		uint64 y = mi.get_CRT(x[i]);

		// Xm += y / mi  mod m
		uint64 Xim = m.set(y);
		uint64 mi_inv = m.inv(m.set(mi.modulus()));
		Xm = m.add(Xm, m.mul(Xim, mi_inv));

		// X += y / mi  mod 1
		X += static_cast<double>(y) / mi.modulus();
		if(X >= 1.0)
		{
			wraps++;
			X -= 1.0;
		}
	}
	cout << " W = " << wraps << endl;

	// Xm -= wraps mod m
	Xm = m.sub(Xm, m.set(wraps));

	// Xm *= M mod m
	Xm = m.mul(Xm, Mm);

	return Xm;
}
