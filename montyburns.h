#pragma once
#include <monty.h>
#include <burns.h>
#include <gmp.h>

class montyburns_cache
{
	public:
		const burns& burns_for_bits(int bits);
		const monty& field(int i);
		const burns& rns(int i);
		
	private:
		void create_montys(int n);
		void create_burns(int n);
		vector<monty> fields;
		vector<burns> rnss;
};

extern montyburns_cache mb;

inline const monty& montyburns_cache::field(int i)
{
	return fields[i];
}

inline const burns& montyburns_cache::rns(int i)
{
	return rnss[i];
}

