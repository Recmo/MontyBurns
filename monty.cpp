#include <monty.h>
#include <primes.h>
#include <iostream>

using namespace std;

monty::monty(uint64 modulus)
{
	assert(modulus > (1ul << 63));
	assert(is_prime(modulus));
	m = modulus;
	
	// φ(2⁶⁴) = 2⁶³
	// k = R^(2⁶³-1)  (implicit mod 2⁶⁴)
	k = 1;
	uint64 b = monty_R();
	for(uint64 e = (1ul << 63) - 1; e; e >>= 1)
	{
		if(e & 1) k *= b;
		b *= b;
	}
	
	two_inv = inv(set_small(2));
	
#ifdef discrete_logarithm
	// Factor φ(m) = m - 1
	phi_factors = prime_factors(m - 1);
	
	// Find the smalest generator
	for(uint64 i = set_small(2); true; i = add(i, one()))
	{
		if(order(i) == m -1)
		{
			g[0] = i;
			break;
		}
	}
	
	// Calculate the powers of the generator
	for(int i = 1; i < 64; i++)
	{
		g[i] = mul(g[i - 1], g[i - 1]);
	}
#endif
	
	/*
	cout << " m = " << m << endl;
	cout << " R = " << R << endl;
	cout << " r = " << r << endl;
	cout << " k = " << k << endl;
	cout << " F = " << phi_factors << endl;
	cout << " g = [";
	for(int i=0; i < 64; i++)
	{
		cout << g[i];
		if(i < 63) cout << ", ";
	}
	cout << "]" << endl;
	cout << endl;
	*/
}

/// Modular exponentiation
/// @param b The base factor in Monty form
/// @param e The exponent
/// @returns The power in Monty form: b^e 2⁻⁶⁴ mod m
uint64 monty::pow(uint64 b, uint64 e) const
{
	// TODO: pow might benefit from inlining
	// constant e allows the compiler to unroll
	// the loop
	
	uint64 p = one();
	for(; e; e >>= 1)
	{
		if(e & 1)
		{
			p = mul(p, b);
		}
		b = mul(b, b);
	}
	return p;
}

#ifdef discrete_logarithm

uint64 monty::order(uint64 n) const
{
	// Start with order o = φ(m) = m - 1
	uint64 o = m - 1;
	
	// The order must be a divisor of φ(m)
	// so try to remove factors of φ(m)
	for(int i=0; i < phi_factors.size(); i++)
	{
		uint64 f = phi_factors[i];
		while(o % f == 0)
		{
			// See if we can remove this factor
			uint64 try_o = o / f;
			if(pow(n, try_o) == one())
			{
				o = try_o;
			}
			else
			{
				break;
			}
		}
	}
	return o;
}

/// Exponentiation of the smallest generator
uint64 monty::exp(uint64 e) const
{
	uint64 p = one();
	for(int i = 0; e; e >>= 1, ++i)
	{
		if(e & 1)
		{
			p = mul(p, g[i]);
		}
	}
	return p;
}

/// Discrete logarithm over the smallest generator
/// @returns e such that g^e = a mod m
uint64 monty::log(uint64 a) const
{
	return log(a, generator());
}

/// General discrete logarithm modulo m
/// @returns e such that a^e = b mod m
uint64 monty::log(uint64 a, uint64 b) const
{
	/// TODO: implement
	return 0;
}

#endif
