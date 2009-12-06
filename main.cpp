#include <stdarg.h>
#include <burns.h>
#include <monty.h>
#include <primes.h>
#include <stdlib.h>
#include <time.h>

double timems(const timespec& start, const timespec& finish)
{
	return (finish.tv_sec - start.tv_sec) * 1000.0 +
	(finish.tv_nsec - start.tv_nsec) / 1e6;
}

void benchmark(int bits)
{
	cout << "{" << bits << ", ";

	// Timers
	clockid_t cpuclock;
	clock_getcpuclockid(0, &cpuclock);
	timespec start, finish;

	// Initialize a Monty Burns system
	burns   mb(bits);

	// Initialize GMP random state
	gmp_randstate_t rnd;
	gmp_randinit_default(rnd);

	// Initialize GMP variables
	mpz_t a, b, n;
	mpz_init(a);
	mpz_init(b);
	mpz_init(n);
	mpz_urandomb(a, rnd, 1 * bits / 2);
	mpz_urandomb(b, rnd, 1 * bits / 2);

	// Calculate the product using GMP
	clock_gettime(cpuclock, &start);
	mpz_mul(n, a, b);
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << ", ";

	// Initialize the MB variables
	vector<uint64> a_mb = mb.set(a);
	vector<uint64> b_mb = mb.set(b);
	vector<uint64> n_mb(mb.size());

	clock_gettime(cpuclock, &start);
	for(int i=0; i < mb.size(); i++)
	{
		const monty& m = mb.montys()[i];
		n_mb[i] = m.mul(a_mb[i], b_mb[i]);
	}
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << "}" << endl;

	// Compare results
	mpz_t n2;
	mpz_init(n2);
	mb.get(n2, n_mb);

	if(mpz_cmp(n, n2) != 0)
	{
		cout << "DIFFERENCE!" << endl;
	}
}

int main()
{
	cout << "{" << endl;
	for(int b=1; b < 10000; b++)
	{
		if (b != 1) cout << ",";
		benchmark(b);
	}
	cout << "}" << endl;

	return 0;
}
