#include <stdarg.h>
#include <burns.h>
#include <monty.h>
#include <primes.h>
#include <stdlib.h>

uint64 factorial(const monty m, uint64 a)
{
	uint64 rm = m.set(1);
	uint64 am = m.set(a);
	uint64 minone = m.neg(m.set(1));
	while(a--)
	{
		rm = m.mul(rm, am);
		am = m.add(am, minone);
	}
	return rm;
}

int main()
{
	// Construct a RNS for a million decimals
	//burns b(4000000);
	burns   b(180);
	int k = 40;
	gmp_printf("M = %Zd\n", b.modulus());

	cout << "Calculating..." << endl;
	vector<uint64> residues;
	residues.reserve(b.size());
	for(int i=0; i < b.size(); i++)
	{
		monty m = b.montys()[i];
		uint64 residue = factorial(m, k);

		// uint64 pow2 = m.pow(m.set(2), 502);
		// pow2 = m.inv(pow2);
		// residue = m.mul(residue, pow2);

		residues.push_back(residue);
		if(i % 1000 == 0) cout << i << " of " << b.size() << endl;
	}

	cout << endl << endl;

	cout << "Converting..." << flush;
	mpz_t n;
	mpz_init(n);
	b.get(n, residues);
	cout << endl;
	gmp_printf("X = %Zd\n", n);

	monty m(9223372036854775837ul);

	double Inv64 = 5.421010862427522E-20;

	cout.precision(16);
	cout << " W   = " << b.count_wraps(residues) << endl;
	//cout << " X/M = " << b.fractional64(residues) << endl;
	//cout << " X/M = " << b.fractional64(residues) << endl;
	cout << " X/M = " << b.fractional(residues) * Inv64 << endl;
	//cout << " W = " << b.count_wraps(residues) << endl;
	//cout << " X mod " << m.modulus() << " = " << m.get(b.mod(residues, m)) << endl;
	cout << " X mod 2⁶⁴ = " << b.mod64(residues) << endl;

	cout << endl << endl;

	return 0;
}
