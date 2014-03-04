#include <montyburns.h>
#include <primes.h>
#include <gmp.h>

montyburns_cache mb;

const burns& montyburns_cache::burns_for_bits(int bits)
{
	// The number of bits gained per prime is
	// more than 63.999999999999 for the first
	// hundred thousand primes.
	// Take one bit extra to accomodate for this
	// slight deviation from 64.
	uint64 n = (bits + 1) / 64 + 1;
	if(n >= 2) n--;
	
	// Ensure sufficient montys are available
	create_montys(n);
	
	if(rnss.size() < n)
	{
		rnss.resize(n);
	}
	if(rnss[n - 1].size() == 0)
	{
		rnss[n - 1] = burns(n);
	}
	
	// Return the burns
	return rnss[n - 1];
}

void montyburns_cache::create_montys(int n)
{
	if(fields.size() >= n) return;
	
	// Find the primes and construct the Montgomery fields
	fields.reserve(n);
	uint64 p = max_uint64;
	if(fields.size())
	{
		p = fields.back().modulus() - 2;
	}
	for(; fields.size() < n; p -= 2)
	{
		if(is_prime(p))
		{
			fields.push_back(monty(p));
		}
	}
}

void montyburns_cache::create_burns(int n)
{
	if(rnss.size() >= n) return;
	
	// Find the primes and construct the Montgomery fields
	rnss.reserve(n);
	for(int i = rnss.size() + 1; i <= n; i++)
	{
		rnss.push_back(burns(i));
	}
}
