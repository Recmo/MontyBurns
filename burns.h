#include <utilities.h>
#include <gmp.h>

class monty;

// Residue Number System
class burns
{
	public:
		burns(int bits);

		const int size() const;
		const vector<monty>& montys() const;

		vector<uint64> set(mpz_t integer) const;
		void get(mpz_t integer, vector<uint64> residues) const;

		uint64 mod(vector<uint64> residues, const monty& m) const;
		uint64 mod64(vector<uint64> residues) const;

		// Returns X mod 10^19
		uint64 mod10(vector<uint64> residues) const;
		uint64 to_uint64(vector<uint64> residues) const;

		const mpz_t& modulus() { return M; }

		uint64 fractional(vector<uint64> residues) const;

		uint64 count_wraps(vector<uint64> residues) const;

	private:
		vector<monty> fields;
		vector<uint64> mrc;
		uint64 Mmod3;

		mpz_t M;
};

inline const int burns::size() const
{
	return fields.size();
}

inline const vector<monty>& burns::montys() const
{
	return fields;
}
