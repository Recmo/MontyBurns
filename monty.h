#pragma once
#include <utilities.h>

// Montgomery Reduced Galois Field
class monty
{
	public:
		monty(uint64 modulus);
		
		// Identity elements in the field
		uint64 zero() const;
		uint64 one() const;
		uint64 two() const;
		
		// Convert to/from uint64 values
		uint64 set(uint64 a) const;
		uint64 get(uint64 a) const;
		
		uint64 add(uint64 a, uint64 b) const;
		uint64 neg(uint64 a) const;
		uint64 sub(uint64 a, uint64 b) const;
		
		uint64 mul(uint64 a, uint64 b) const;
		uint64 pow(uint64 a, uint64 b) const;
		uint64 inv(uint64 a) const;
		
		uint64 order(uint64 n) const;
		
		uint64 modulus() const;
		uint64 monty_R() const;
		uint64 monty_r() const;
		uint64 monty_k() const;
		uint64 genertor() const;
		
		// Non Montgomery multiplication and power
		uint64 modular_mul(uint64 a, uint64 b) const;
		
	private:
		uint64 m; // The modulus
		uint64 R; // = 2⁶⁴ mod m = 2⁶⁴ - m
		uint64 r; // = 2⁻⁶⁴ mod m = R^(m-2) mod m
		uint64 k; // = (-m)⁻¹ mod 2⁶⁴ = R^(2⁶⁴-2) mod 2⁶⁴
		vector<uint64> phi_factors; // Prime factors of m - 1
		uint64 g; // Smallest generator
};

/// Converts from integer to Montgomery reduced form
/// @param a The number in a mod m form
/// @returns The Monty form: a 2⁶⁴ mod m
inline uint64 monty::set(uint64 a) const
{
	if(fast_false(a >= m)) a -= m;
	return modular_mul(a, R);
}

/// Converts from Montgomery reduced form to modular
/// @param a The number Monty form: a 2⁶⁴ mod m  form
/// @returns The modular form: a mod m
inline uint64 monty::get(uint64 a) const
{
	return modular_mul(a, r);
}

inline uint64 monty::zero() const
{
	return 0;
}

inline uint64 monty::one() const
{
	return R;
}

inline uint64 monty::two() const
{
	// R << m   so no reduction is necesary
	return R + R;
}


/// Modular addition
/// @param a The first term, a < m
/// @param a The second term, b < m
/// @return The sum, a b mod m
inline uint64 monty::add(uint64 a, uint64 b) const
{
	uint64 s = a + b;
	if(s < a) s += R;
	if(fast_false(s >= m)) s -= m;
	return s;
}

/// Modular negation
/// @param a The number, a < m
/// @return The sum, a b mod m
inline uint64 monty::neg(uint64 a) const
{
	return m - a;
}

/// Modular substraction
/// @param a The first term, a < m
/// @param a The second term, b < m
/// @return The substraction, a - b mod m
inline uint64 monty::sub(uint64 a, uint64 b) const
{
	return (a > b) ? a - b : (m - b) + a;
}

/// Multiplication with Montgomery reduction
/// @param a The first factor in Monty form
/// @param b The second factor in Monty form
/// @returns The product in Monty form: a b 2⁻⁶⁴ mod 2⁶⁴
inline uint64 monty::mul(uint64 a, uint64 b) const
{
	// h 2⁶⁴ + l = a b
	uint64 h, l;
	asm("mulq %3" : "=a"(l), "=d"(h) : "a"(a), "rm"(b));
	
	// If it is a multiple of 2⁶⁴
	// then return the quotient
	if(fast_false(l == 0)) return h;
	
	// u = t n  mod  r
	//   = (r mod r) n mod r
	//   = (l1 | l2) n mod r
	uint64 u = l * k;
	
	// vh 2⁶⁴ + vl = u m
	uint64 vh, vl;
	asm("mulq %3;" : "=a"(vl), "=d"(vh) : "a"(u), "rm"(m));
	
	// a = vh + h + 1 (take care of overflow)
	uint64 p = h + vh + 1;
	if(p < h) p += R;
	if(fast_false(p >= m)) p -= m;
	
	return p;
}

/// Modular exponentiation
/// @param a The number to invert in Monty form
/// @returns The multiplicative inverse in Monty form: a⁻^1 2⁻⁶⁴ mod m
inline uint64 monty::inv(uint64 a) const
{
	// Apply Fermat's little theorem
	// φ(m) = m - 1  (m is prime)
	return pow(a, m - 2);
}


inline uint64 monty::modulus() const
{
	return m;
}

inline uint64 monty::monty_R() const
{
	return R;
}

inline uint64 monty::monty_r() const
{
	return r;
}

inline uint64 monty::monty_k() const
{
	return k;
}

/// Modular multiplication
/// @param a The first factor, a < m
/// @param a The second factor, b < m
/// @param m The modulus
/// @return The reduced product, a b mod m < m
inline uint64 monty::modular_mul(uint64 a, uint64 b) const
{
	// Perform 128 multiplication and division
	uint64 q; // q = ⌊a b / m⌋
	uint64 r; // r = a b mod m
	asm(
	"mulq %3;"
	"divq %4;"
	: "=a"(q), "=d"(r)
	: "a"(a), "rm"(b), "rm"(m));
	return r;
}
