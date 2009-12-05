#include <primes.h>
#include <stdlib.h>

const uint64s small_primes = prime_sieve(1<<20);

const prime_wheel wheel;

uint64s prime_sieve(uint64 max)
{
	vector<bool> sieve((max + 1) / 2);
	uint64s primes;
	primes.push_back(2);
	for(uint64 i = 3; i <= max; i += 2)
	{
		if(sieve[(i - 3) / 2]) continue;
		for(uint64 m = i * i; m <= max; m += 2 * i)
			sieve[(m - 3) / 2] = true;
		primes.push_back(i);
	}
	return primes;
}

prime_wheel::prime_wheel()
{
	uint64 max_size = small_primes[small_primes.size() - 1];

	max_size = 211;

	size = 1;
	for(factors = 0; size < max_size; factors++)
	{
		size *= small_primes[factors];
	}
	factors--;
	size /= small_primes[factors];

	uint64 previous = 1;
	for(int i = 0; i < size; i++)
	{
		bool include = true;
		for(int j = 0; j < factors; j++)
		{
			if((i % small_primes[j] == 0))
			{
				include = false;
			}
		}
		if(!include) continue;
		deltas.push_back(i - previous);
		previous = i;
	}
	deltas.push_back(size - previous + 1);
}

vector<uint64> prime_factors_wheel(uint64 n)
{
	vector<uint64> f;

	// Remove the twos
	if((n & 1) == 0)
	{
		f.push_back(2);
		do n >>= 1; while((n & 1) == 0);
	}

	// Remove prime wheel factors
	uint64 p = 2;
	for(int i = 1; p < wheel.size; i++)
	{
		p = small_primes[i];
		if(fast_false(n % p == 0))
		{
			f.push_back(p);
			do n /= p; while(n % p == 0);
		}
	}

	// Turn the prime wheel
	for(p = wheel.size + 1; p * p < n;)
	{
		for(int j = 0; j < wheel.deltas.size(); j++)
		{
			if(fast_false(n % p == 0))
			{
				f.push_back(p);
				do n /= p; while(n % p == 0);
			}
			p += wheel.deltas[j];
		}
	}

	// What's left must be prime or one
	if(n != 1) f.push_back(n);

	return f;
}

uint64 mul(uint64 a, uint64 b, uint64 m)
{
	// Perform 128 multiplication and division
	uint64 q; // q = ⌊a b / m⌋
	uint64 r; // r = a b mod m
	asm(
		"mulq %3;"
		"divq %4"
		: "=a"(q), "=d"(r)
		: "a"(a), "rm"(b), "rm"(m));
		return r;
}

uint64 pow(uint64 b, uint64 e, uint64 m)
{
	uint64 r = 1;
	for(; e; e >>= 1)
	{
		if(e & 1) r = mul(r, b, m);
		b = mul(b, b, m);
	}
	return r;
}

/// Miller-Rabin probabilistic primality testing
/// @param n The number to test for  primality
/// @param k The witness for primality
/// @returns True iff when n is a k-stong pseudoprime
bool miller_rabin(uint64 n, uint64 k)
{
	// Factor n-1 as d*2^s
	uint64 s = 0;
	uint64 d = n - 1;
	for(; !(d & 1); s++)
		d >>= 1;

	// Verify x = k^(d 2^i) mod n != 1
	uint64 x = pow(k % n, d, n);
	if(x == 1 || x == n-1) return true;
	while(s-- > 1)
	{
		// x = x^2 mod n
		x = mul(x, x, n);
		if(x == 1) return false;
		if(x == n-1) return true;
	}
	return false;
}

/// Miller-Rabin probabilistic primality testing
/// @param n The number to test for  primality
/// @returns False when n is not a prime
bool is_prime(uint64 n)
{
	// Handle small primes fast
	for(int i = 0; i < 100; i++)
	{
		uint64 p = small_primes[i];
		if(n == p) return true;
		if(n % p == 0) return false;
	}

	// Do a few Miller-Rabin rounds
	for(int i = 0; i < 10; i++)
	{
		if (!miller_rabin(n, small_primes[i])) return false;
	}
	return true;
}

// Source: http://en.wikipedia.org/wiki/Binary_GCD_algorithm
uint64 gcd(uint64 u, uint64 v)
{
	int shift;

	// GCD(0,x) := x
	if (u == 0 || v == 0)
		return u | v;

	/* Let shift := lg K, where K is the greatest power of 2
	dividing both u and v. */
	for (shift = 0; ((u | v) & 1) == 0; ++shift)
	{
		u >>= 1;
		v >>= 1;
	}

	while ((u & 1) == 0)
		u >>= 1;

	/* From here on, u is always odd. */
	do
	{
		while ((v & 1) == 0)  /* Loop X */
			v >>= 1;

		/* Now u and v are both odd, so diff(u, v) is even.
		Let u = min(u, v), v = diff(u, v)/2. */
		if (u < v)
		{
			v -= u;
		}
		else
		{
			uint64 diff = u - v;
			u = v;
			v = diff;
		}
		v >>= 1;
	} while (v != 0);

	return u << shift;
}

void brents_factor(uint64 n, vector<uint64>& f)
{
	if(is_prime(n))
	{
		f.push_back(n);
		return;
	}

	uint64 a = random() % n;
	uint64 x = random() % n;
	uint64 y = 1;

	for(int i = 0; i*i / 2 < n; i++)
	{
		// x = x² + a mod n
		x = mul(x, x, n);
		x += a;
		if(x < a) x += (max_uint64 - n) + 1;
		x %= n;

		uint64 g = gcd(n, y - x);
		if(g != 1 && g != n)
		{
			n /= g;
			brents_factor(g, f);
			if(n != g) brents_factor(n, f);
			return;
		}

		if((i & (i-1)) == 0) y = x;
	}

	// Found no factors, yet n is not a prime, retry
	brents_factor(n, f);
}

vector<uint64> prime_factors(uint64 n)
{
	vector<uint64> f;

	// Remove factors of two
	if((n & 1) == 0)
	{
		f.push_back(2);
		do n >>= 1; while((n & 1) == 0);
	}
	if(n == 0 || n == 1) return f;

	// Find small factors
	for(int i=1; i < 100; i++)
	{
		uint64 p = small_primes[i];
		if(n % p == 0)
		{
			f.push_back(p);
			do n /= p; while(n % p ==0);
		}
	}
	if(n == 1) return f;

	// Call on the recursive Brent's factorization
	brents_factor(n, f);

	// Todo: Sort and remove duplicates

	return f;
}

vector<uint64> prime_factors_old(uint64 n)
{
	vector<uint64> f;

	// Remove the twos
	if((n & 1) == 0)
	{
		f.push_back(2);
		do n >>= 1; while((n & 1) == 0);
	}

	// Remove other factors
	uint64 i;
	for(i=3; i*i < n; i += 2)
	{
		if(fast_false(n % i == 0))
		{
			f.push_back(i);
			do n /= i; while(n % i == 0);
		}
	}
	f.push_back(n);
	return f;
}
