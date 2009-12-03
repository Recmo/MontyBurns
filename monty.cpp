#include <monty.h>
#include <primes.h>
#include <iostream>

using namespace std;

monty::monty(uint64 modulus)
{
	m = modulus;
	R = max_uint64 - m + 1;
	
	// φ(m) = m - 1  (m is prime)
	// r = R^(m-2) mod m
	r = 1;
	uint64 b = R;
	for(uint64 e = m - 2; e; e >>= 1)
	{
		if(e & 1) r = modular_mul(r, b);
		b = modular_mul(b, b);
	}
	
	// φ(2⁶⁴) = 2⁶³
	// k = R^(2⁶³-1)  (implicit mod 2⁶⁴)
	k = 1;
	b = R;
	for(uint64 e = (1ul << 63) - 1; e; e >>= 1)
	{
		if(e & 1) k *= b;
		b *= b;
	}
	
	/*
	// Factor φ(m) = m - 1
	phi_factors = prime_factors(m - 1);
	
	// Find the smalest generator
	for(uint64 i = two(); true; i = add(i, one()))
	{
		if(order(i) == m -1)
		{
			g = i;
			break;
		}
	}
	*/
	
	/*
	cout << " m = " << m << endl;
	cout << " R = " << R << endl;
	cout << " r = " << r << endl;
	cout << " k = " << k << endl;
	cout << " F = ";
	for(int i=0; i < phi_factors.size(); i++)
	{
		cout << phi_factors[i] << ", ";
	}
	cout << endl;
	cout << " g = " << get(g) << endl;
	cout << endl;
	*/
}

/// Modular exponentiation
/// @param b The base factor in Monty form
/// @param e The exponent
/// @returns The power in Monty form: b^e 2⁻⁶⁴ mod m
uint64 monty::pow(uint64 b, uint64 e) const
{
	// e %= m   (not useful, very unlikely)
	uint64 p = R;
	for(; e; e >>= 1)
	{
		if(e & 1) p = mul(p, b);
		b = mul(b, b);
	}
	return p;
}

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



