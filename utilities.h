#pragma once
#include <iostream>
#include <vector>

using namespace std;

typedef unsigned int uint32;
typedef unsigned long int uint64;
typedef vector<uint64> uint64s;
typedef uint64s::iterator uint64s_i;
typedef uint64s::const_iterator uint64s_ci;

const uint64 max_uint64 = 0xFFFFFFFFFFFFFFFF;

#define fast_true(a) __builtin_expect(a, 1)
#define fast_false(a) __builtin_expect(a, 0)
