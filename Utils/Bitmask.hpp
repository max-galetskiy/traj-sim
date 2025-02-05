#pragma once

/*
    Simple, light-weight bitmask implementation to use instead of built-in lookup-containers like std::unordered_map or boolean arrays
 */
class Bitmask {

    private:

        int nr_bitmasks = 1;

        static const int bitmask_size = 8 * sizeof(unsigned long long);

        unsigned long long* bitmasks;

    public:

        Bitmask(int required_bits);
        ~Bitmask(){delete[] bitmasks;}

        bool bit_set(int b);
        void set_bit(int b);

};
