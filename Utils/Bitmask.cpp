#include <cmath>
#include "Bitmask.hpp"

Bitmask::Bitmask(int required_bits) {
    nr_bitmasks = std::ceil(required_bits * 1.f / bitmask_size);
    bitmasks = new unsigned long long[nr_bitmasks];

    for(int i = 0; i < nr_bitmasks; i++){
        bitmasks[i] = 0;
    }
}

bool Bitmask::bit_set(int b) {

    if(b < bitmask_size){
        return bitmasks[0] & (1ULL << b);
    }

    return bitmasks[b / bitmask_size] & (1ULL << (b % bitmask_size));

}

void Bitmask::set_bit(int b) {

    if(b < bitmask_size){
        bitmasks[0] |= (1ULL << b);
        return;
    }

    bitmasks[b / bitmask_size] |= (1ULL << (b % bitmask_size));

}
