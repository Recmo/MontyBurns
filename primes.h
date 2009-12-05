#include <utilities.h>

/// The ubiquitous table of small primes
extern const uint64s small_primes;

/// Erastotheles Prime Sieve
/// @param max The largest number to consider
/// @return A complete list of primes smaller than max
uint64s prime_sieve(uint64 max);

uint64 pollard_factor(uint64 n);

vector<uint64> prime_factors(uint64 n);
vector<uint64> factors(uint64 n);

/// Miller-Rabin probabilistic primality testing
/// @param n The number to test for  primality
/// @param k The witness for primality
/// @returns True iff when n is a k-stong pseudoprime
bool miller_rabin(uint64 n, uint64 k);

/// Miller-Rabin probabilistic primality testing
/// @param n The number to test for  primality
/// @returns False when n is not a prime
bool is_prime(uint64 n);

uint64 gcd(uint64 a, uint64 b);


