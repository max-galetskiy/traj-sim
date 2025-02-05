#include "Utilities.hpp"
#include <random>



void inside_out_shuffle(int* arr, unsigned int n){

    // seeding random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    for(int i = 0; i < n; i++){

        std::uniform_int_distribution<> distr(0,i);
        int j = distr(gen);

        if(j != i){
            arr[i] = arr[j];
        }

        arr[j] = i;

    }

}