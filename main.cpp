#include <montyburns.h>
#include <stdarg.h>
#include <burns.h>
#include <monty.h>
#include <primes.h>
#include <stdlib.h>
#include <time.h>

gmp_randstate_t rnd;

double timems(const timespec& start, const timespec& finish)
{
	return (finish.tv_sec - start.tv_sec) * 1000.0 +
	(finish.tv_nsec - start.tv_nsec) / 1e6;
}

void benchmark_init(int bits)
{
	cout << "{" << bits << ", " << flush;
	
	// Timers
	clockid_t cpuclock;
	clock_getcpuclockid(0, &cpuclock);
	timespec start, finish;
	
	// Initialize a Monty Burns system
	clock_gettime(cpuclock, &start);
	const burns& mbrns = mb.burns_for_bits(bits);
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << flush;
	
	cout << "}" << endl;
}

void benchmark_conversions(int bits)
{
	cout << "{" << bits << ", " << flush;
	
	// Timers
	clockid_t cpuclock;
	clock_getcpuclockid(0, &cpuclock);
	timespec start, finish;
	
	// Initialize a Monty Burns system
	const burns& mbrns = mb.burns_for_bits(bits);
	
	// Initialize GMP variables
	mpz_t a, b, n, n2;
	mpz_init(n);
	mpz_init(n2);
	mpz_urandomb(n, rnd, bits);
	
	// Convert to montyburns
	clock_gettime(cpuclock, &start);
	vector<uint64> n_mb = mbrns.set(n);
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << ", " << flush;
	
	// Convert to GMP
	clock_gettime(cpuclock, &start);
	mbrns.get(n2, n_mb);
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << ", " << flush;
	
	// Convert from rns to mixed radix
	clock_gettime(cpuclock, &start);
	mbrns.to_mixed_radix(n_mb);
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << ", " << flush;
	
	// Convert from mixed radix to rns
	clock_gettime(cpuclock, &start);
	mbrns.from_mixed_radix(n_mb);
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << flush;
	
	// Compare results
	if(mpz_cmp(n, n2) != 0)
	{
		cout << "DIFFERENCE!" << endl;
	}
	
	cout << "}" << endl;
}

void benchmark_utils(int bits)
{
	cout << "{" << bits << ", " << flush;
	
	// Timers
	clockid_t cpuclock;
	clock_getcpuclockid(0, &cpuclock);
	timespec start, finish;
	
	// Initialize a Monty Burns system
	const burns& mbrns = mb.burns_for_bits(bits);
	
	// Initialize GMP variables
	mpz_t a, b, n;
	mpz_init(n);
	mpz_urandomb(n, rnd, bits);
	
	// Convert to montyburns
	vector<uint64> n_mb = mbrns.set(n);
	
	// Calculate fractional
	clock_gettime(cpuclock, &start);
	uint64 fract = mbrns.fractional(n_mb);
	//cout << "X/M = (" << fract << " ... " << fract + n_mb.size() << ")/2^64" << endl;
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << ", " << flush;
	
	// Calculate fractional using GMP
	clock_gettime(cpuclock, &start);
	mpz_mul_2exp(n, n, 64);
	mpz_fdiv_q(n, n, mbrns.modulus());
	//gmp_printf("N = %Zd\n", n);
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << ", " << flush;
	
	// Calculate mod64
	clock_gettime(cpuclock, &start);
	uint64 mod64 = mbrns.mod64(n_mb);
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << ", " << flush;
	
	// Calculate mod m
	/*
	clock_gettime(cpuclock, &start);
	mbrns.mod(n_mb, 1233);
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << flush;
	*/
	
	cout << "}" << endl;
	
	cout << "mod64 = " << mod64 << endl;
	cout << "mod64 = " << n->_mp_d[0] << endl;
}

void benchmark_mul(int bits)
{
	cout << "{" << bits << ", ";

	// Timers
	clockid_t cpuclock;
	clock_getcpuclockid(0, &cpuclock);
	timespec start, finish;

	// Initialize a Monty Burns system
	const burns& mbrns = mb.burns_for_bits(bits);

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
	vector<uint64> a_mb = mbrns.set(a);
	vector<uint64> b_mb = mbrns.set(b);
	vector<uint64> n_mb(mbrns.size());

	clock_gettime(cpuclock, &start);
	for(int i=0; i < mbrns.size(); i++)
	{
		const monty& m = mb.field(i);
		n_mb[i] = m.mul(a_mb[i], b_mb[i]);
	}
	clock_gettime(cpuclock, &finish);
	cout << timems(start, finish) << "}" << endl;

	// Compare results
	mpz_t n2;
	mpz_init(n2);
	mbrns.get(n2, n_mb);

	if(mpz_cmp(n, n2) != 0)
	{
		cout << "DIFFERENCE!" << endl;
	}
}

int main()
{
	// Initialize GMP random state
	gmp_randinit_default(rnd);
	
	cout << "{" << endl;
	for(int b=1; b < 1<<21; b <<= 1)
	{
		if (b != 1) cout << ",";
		benchmark_init(b);
	}
	cout << "}" << endl;
	
	/*
	cout << "{" << endl;
	for(int b=1; b < 1<<21; b <<= 1)
	{
		if (b != 1) cout << ",";
		benchmark_init(b);
	}
	cout << "}" << endl;
	
	cout << "{" << endl;
	for(int b=1; b < 1<<21; b <<= 1)
	{
		if (b != 1) cout << ",";
		benchmark_conversions(b);
	}
	cout << "}" << endl;
	*/
	
	cout << "{" << endl;
	for(int b=1; b < 1<<21; b <<= 1)
	{
		if (b != 1) cout << ",";
		benchmark_utils(b);
	}
	cout << "}" << endl;
	
	
	
	return 0;
}
