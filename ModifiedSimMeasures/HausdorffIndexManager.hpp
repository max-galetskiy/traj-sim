#pragma once
#include <vector>
#include <random>

#include "../Utils/Bitmask.hpp"

class HausdorffIndexManager{
    public:

        void insert_blacklisted(int i, int j){
            if (!(blacklisted_i_values.bit_set(i))) {nr_blacklisted_i++;}
            if (!(blacklisted_j_values.bit_set(j))){ nr_blacklisted_j++;}
            blacklisted_i_values.set_bit(i);
            blacklisted_j_values.set_bit(j);
        }

        void insert_prioritized(int i, bool i_value){   // because only entirely grey rows/columns are prioritized, it could happen that row i is entirely grey, while column j isnt --> therefore only i or j are inserted and not both

            (i_value) ? prioritized_i_values.emplace_back(i) : prioritized_j_values.emplace_back(i);
            (i_value) ? (is_prioritized_i.set_bit(i)) : (is_prioritized_j.set_bit(i));
        }

        void insert_unprioritized(int i, int j){        // because we unprioritize rows/columns that contain a white point, both row i and column j will be unprioritized when finding the point (i,j)

            if(!is_unprioritized_i.bit_set(i)){
                unprioritized_i_values.emplace_back(i);
                is_unprioritized_i.set_bit(i);
            }
            if(!is_unprioritized_j.bit_set(j)){
                unprioritized_j_values.emplace_back(j);
                is_unprioritized_j.set_bit(j);
            }

        }

        void toggle_type(bool consider_i_values){
            consider_i = consider_i_values;
        }

        void decide_hausdorff(bool lt_epsilon){
            hausdorff_decided = true;
            hausdorff_LT_epsilon = lt_epsilon;
        }

        bool in_blacklisted(int x){
            return (consider_i) ? (blacklisted_i_values.bit_set(x)) : (blacklisted_j_values.bit_set(x));
        }

        bool in_prioritized(int x){
            return (consider_i) ? (is_prioritized_i.bit_set(x)) : (is_prioritized_j.bit_set(x));
        }

        bool in_unprioritized(int x){
            return (consider_i) ? (is_unprioritized_i.bit_set(x)) : (is_unprioritized_j.bit_set(x));
        }

        bool decided() {return hausdorff_decided;}
        bool lt_eps() { return hausdorff_LT_epsilon;}

        std::vector<int> * get_prioritized(){
            return (consider_i ? &prioritized_i_values : &prioritized_j_values);
        }

        std::vector<int> * get_unprioritized(){
            return (consider_i ? &unprioritized_i_values : &unprioritized_j_values);
        }

        HausdorffIndexManager(unsigned int n, unsigned int m) : blacklisted_i_values(Bitmask(n)), blacklisted_j_values(Bitmask(m)),
                                                                is_prioritized_i(Bitmask(n)), is_prioritized_j(Bitmask(m)),
                                                                is_unprioritized_i(Bitmask(n)), is_unprioritized_j(Bitmask(m))
        {
            prioritized_i_values.reserve(n);
            prioritized_j_values.reserve(m);

            unprioritized_i_values.reserve(n);
            unprioritized_j_values.reserve(m);
        }

        int nr_blacklisted_i = 0;
        int nr_blacklisted_j = 0;

    private:

        bool consider_i = true;

        bool hausdorff_decided = false;
        bool hausdorff_LT_epsilon;



        Bitmask blacklisted_i_values;
        Bitmask blacklisted_j_values;

        std::vector<int> prioritized_i_values;
        std::vector<int> prioritized_j_values;

        Bitmask is_prioritized_i;
        Bitmask is_prioritized_j;

        std::vector<int> unprioritized_i_values;
        std::vector<int> unprioritized_j_values;

        Bitmask is_unprioritized_i;
        Bitmask is_unprioritized_j;

};

void shuffle_index(int* arr, int i, unsigned int size, std::mt19937& gen){
    std::uniform_int_distribution<> distr(0,size);
    int j = distr(gen);

    if(j != i){
        arr[size] = arr[j];
    }

    arr[j] = i;
}

int modified_inside_out_shuffle(int* arr, int* prior, int* unprior, unsigned int n, HausdorffIndexManager& hd){

    // seed rng
    std::random_device rd;
    std::mt19937 gen(rd());

    int size = 0;
    int prior_size = 0;
    int unprior_size = 0;

    for(int i = 0; i < n; i++){

        if(!(hd.in_blacklisted(i) || hd.in_prioritized(i) || hd.in_unprioritized(i))){
            shuffle_index(arr, i, size, gen);
            size++;
        }
        else if(hd.in_prioritized(i)){
            shuffle_index(prior, i, prior_size, gen);
            prior_size++;
        }
        else if(hd.in_unprioritized(i)){
            shuffle_index(unprior, i, unprior_size, gen);
            unprior_size++;
        }

    }

    return size;

}


