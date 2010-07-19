#pragma once
#include <utilities.h>
#include <gmp.h>

class monty;

// Residue Number System
class burns
{
	public:
		burns();
		burns(int size);
		
		const int size() const;
		const vector<monty>& montys() const;
		
		vector<uint64> set_uint64(const uint64 X) const;
		vector<uint64> set(const uint64* bits, int length) const;
		vector<uint64> set(mpz_t integer) const;
		void get(mpz_t integer, const vector<uint64>& residues) const;
		
		void to_mixed_radix(vector<uint64>& x) const;
		void from_mixed_radix(vector<uint64>& x) const;
		
		bool equals(const vector<uint64>& x, uint64 y) const;
		bool equals(const vector<uint64>& x, const vector<uint64>& y) const;
		int compare(const vector<uint64>& x, const vector<uint64>& y);
		
		
		uint64 mod(const vector<uint64>& residues, const monty& m) const;
		uint64 mod64(const vector<uint64>& residues) const;
		
		// Returns X mod 10^19
		uint64 mod10(const vector<uint64>& residues) const;
		uint64 to_uint64(const vector<uint64>& residues) const;
		
		void safe_div(vector<uint64>& x, const vector<uint64>& y) const;
		void safe_div(vector<uint64>& x, const uint64 y, uint64& remainder) const;
		
		bool can_shrink(const vector<uint64>& residues) const;
		
		const mpz_t& modulus() { return M; }
		
		uint64 fractional(const vector<uint64>& residues) const;
		uint64 count_wraps(const vector<uint64>& residues) const;
		
	private:
		int nsize;
		vector<uint64> mrc; // = (r_i m_i / M) mod m_i (non reduced)
		
		vector<vector<uint64> > mji;
		
		uint64 Mmod3; // = M mod 3
		mpz_t M; // = m_1 * m_2 ... m_n
};

inline const int burns::size() const
{
	return nsize;
}

