#include <primes.h>
#include <stdlib.h>

const uint64s small_primes = prime_sieve(1<<20);

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
	for(; e; e >>= 1) {
		if(e & 1)
			r = mul(r, b, m);
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
	if(x == 1 || x == n-1)
		return true;
	while(s-- > 1) {
		// x = x^2 mod n
		x = mul(x, x, n);
		if(x == 1)
			return false;
		if(x == n - 1)
			return true;
	}
	return false;
}

/// Miller-Rabin probabilistic primality testing
/// @param n The number to test for  primality
/// @returns False when n is not a prime
bool is_prime(uint64 n)
{
	// Handle small primes fast
	for(int i = 0; i < 100; i++) {
		uint64 p = small_primes[i];
		if(n == p)
				return true;
		if(n % p == 0)
			return false;
	}

	// Do a few Miller-Rabin rounds
	for(int i = 0; i < 10; i++) {
		if (!miller_rabin(n, small_primes[i])) return false;
	}
	return true;
}

// Source: http://en.wikipedia.org/wiki/Binary_GCD_algorithm
uint64 gcd(uint64 u, uint64 v)
{
	if (u == 0 || v == 0) return u | v;

	// Remove common twos in u and v
	int shift;
	for (shift = 0; ((u | v) & 1) == 0; ++shift)
	{
		u >>= 1;
		v >>= 1;
	}

	// Remove twos from u
	while ((u & 1) == 0) u >>= 1;

	do
	{
		// Remove twos in v
		while ((v & 1) == 0) v >>= 1;

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

uint64 small_factor(uint64 n)
{
	if((n & 1) == 0) return 2;
	for(uint64 i = 3; ; i += 2)
	{
		if(n % i == 0) return i;
	}
}

uint64 brent_pollard_factor(uint64 n)
{
	const uint64 m = 1000;
	uint64 a, x, y, ys, r, q, g;
	do a = random() % n; while(a==0||a==n-2);
	y = random() % n;
	r = 1;
	q = 1;

	do
	{
		x = y;
		for(uint64 i=0; i < r; i++)
		{
			// y = y² + a mod n
			y = mul(y, y, n);
			y += a;
			if(y < a) y += (max_uint64 - n) + 1;
			y %= n;
		}

		uint64 k =0;
		do
		{
			for(uint64 i=0; i < m && i < r-k; i++)
			{
				ys = y;

				// y = y² + a mod n
				y = mul(y, y, n);
				y += a;
				if(y < a) y += (max_uint64 - n) + 1;
				y %= n;

				// q = q |x-y| mod n
				q = mul(q, (x>y)?x-y:y-x, n);
			}
			g = gcd(q, n);
			k += m;
		} while(k < r && g == 1);

		r <<= 1;
	} while(g == 1);

	if(g == n)
	{
		do
		{
			// ys = ys² + a mod n
			ys = mul(ys, ys, n);
			ys += a;
			if(ys < a) ys += (max_uint64 - n) + 1;
			ys %= n;

			g = gcd((x>ys)?x-ys:ys-x, n);
		} while(g == 1);
	}

	return g;
}

vector<uint64> prime_factors(uint64 n)
{
	vector<uint64> factors;
	vector<uint64> primes;

	uint64 factor = brent_pollard_factor(n);
	factors.push_back(n / factor);
	factors.push_back(factor);

	do
	{
		uint64 m = factors[factors.size() - 1];
		factors.pop_back();

		if(m == 1) continue;

		if(is_prime(m))
		{
			primes.push_back(m);

			// Remove the prime from the other factors
			for(int i=0; i < factors.size(); i++)
			{
				uint64 k = factors[i];
				if(k % m == 0)
				{
					do k /= m; while(k % m == 0);
					factors[i] = k;
				}
			}
		}
		else
		{
			factor = (m < 100) ? small_factor(m) : brent_pollard_factor(m);
			factors.push_back(m / factor);
			factors.push_back(factor);
		}
	} while(factors.size());

	return primes;
}
