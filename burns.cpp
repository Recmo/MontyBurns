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
	
	/*
	for(int i=0; i < fields.size(); i++)
	{
		cout << " m_" << i << " = " << fields[i].modulus() << endl;
	}
	*/
	
	// Construct M and M_approx
	mpz_init(M);
	mpz_set_ui(M, 1);
	for(int i=0; i < n; i++)
	{
		uint64 modulus = fields[i].modulus();
		mpz_mul_ui(M, M, modulus);
	}
	// gmp_printf(" M = %Zd\n", M);
	
	// Construct the M inverses
	Minv.reserve(n);
	for(int i=0; i < n; i++)
	{
		const monty& m = fields[i];
		uint64 M = m.one();
		for(int j = 0; j < n; j++)
		{
			if(i == j) continue;
			uint64 modulus = m.set(fields[j].modulus());
			M = m.mul(M, modulus);
		}
		Minv.push_back(m.inv(M));
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
		
		// y = x_i  Minv_i  mod m_i
		uint64 y = m.get(m.mul(x[i], Minv[i]));
		
		// X += Mi y
		mpz_addmul_ui(X, Mi, y);
	}
	mpz_clear(Mi);
	
	// Reduce mod M
	// mpz_mod(X, X, M);
}

uint64 inv64(uint64 b)
{
	uint64 p = 1;
	for(int i=64; i; i--)
	{
		p *= b;
		b *= b;
	}
	return p;
}

double burns::fractional(vector<uint64> x) const
{
	double X = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];
		double y = m.get(m.mul(x[i], Minv[i]));
		X += y / m.modulus();
		if(X >= 1.0) X -= 1.0;
	}
	return X;
}

uint64 burns::fractional64(vector<uint64> x) const
{
	uint64 X = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];
		uint64 y = m.get(m.mul(x[i], Minv[i]));
		X += y * (inv64(m.modulus()) >> 1);
		if(X >= 1.0) X -= 1.0;
	}
	return X;
}

uint64 burns::count_wraps(vector<uint64> x) const
{
	uint64 wraps = 0;
	double X = 0;
	for(int i=0; i < size(); i++)
	{
		const monty& m = fields[i];
		double y = m.get(m.mul(x[i], Minv[i]));
		X += y / m.modulus();
		if(X >= 1.0)
		{
			wraps++;
			X -= 1.0;
		}
	}
	return wraps;
}

uint64 burns::mod64(vector<uint64> x) const
{
	cout << inv64(5) << endl;
	
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
		uint64 y = mi.get(mi.mul(x[i], Minv[i]));
		
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
		uint64 y = mi.get(mi.mul(x[i], Minv[i]));
		
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
	// Xm -= wraps mod m
	Xm = m.sub(Xm ,m.set(wraps));
	
	// Xm *= M mod m
	Xm = m.mul(Xm, Mm);
	
	return Xm;
}
