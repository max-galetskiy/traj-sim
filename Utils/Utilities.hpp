#pragma once
#include <utility>
#include <cstdlib>

// Defines for quicker formatting of output
#define TAB std::setw(LW) << std::left
#define COUT_TAB std::setw(CW) << std::left

/*
	Inside-out variation of the Fisher-Yard shuffle
*/
void inside_out_shuffle(int* arr, unsigned int n);

// Pair hash used for hash-sets of integer pairs
struct pair_hash
{
    size_t operator()(const std::pair<int,int> &x) const
    {
        return x.first ^ x.second;
    }
};