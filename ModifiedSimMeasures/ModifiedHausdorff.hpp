#pragma once
#include <cmath>
#include <algorithm>
#include "../TrajModels/Trajectory.hpp"
#include "../DistCache/DistCache.hpp"
#include "../Utils/Utilities.hpp"
#include "HausdorffIndexManager.hpp"


template<class T>
bool modified_dirHd_LT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, bool flipped, float epsilon, HausdorffIndexManager& hd_ind, std::shared_ptr<DistCache<T>> cache = NULL) {

    // Missing: Remove Intersection

    const unsigned int n = t1.length();
    const unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return (m == n);
    }

    // Reorder Points
    hd_ind.toggle_type(!flipped); // changes whether we're considering indexes of the original t1 or indexes of the original t2

    int* ind_1 = new int[n];
    int* prior = new int[hd_ind.get_prioritized()->size()];
    int* unprior = new int[hd_ind.get_unprioritized()->size()];
    int* ind_2 = new int[m];
    int outer_loop_size = modified_inside_out_shuffle(ind_1,prior, unprior,n,hd_ind);

    inside_out_shuffle(ind_2,m);

    // Calculate directed Hausdorff
    float hmax = 0;

    // Reorder and Iterate over prioritized if not empty
    for(int i = 0; i < hd_ind.get_prioritized()->size(); i++){
        // NN Search
        float hmin = INFINITY;

        for (unsigned int j = 0; j < m; j++) {

            float d;

            if(cache){
                d = flipped ? (*cache).cacheDistance(ind_2[j],prior[i]) : (*cache).cacheDistance(prior[i],ind_2[j]);
            }
            else{
                d = T::dist(t1.at(prior[i]), t2.at(ind_2[j]));
            }

            if (d < hmax) {
                hmin = INFINITY;   // when breaking the loop we have to reset hmin to preserve the integrity of the "hmin > hmax" check
                break;
            }

            if (d < hmin) {
                hmin = d;
            }

        }

        if (hmin > hmax && !(std::isinf(hmin))) {
            hmax = hmin;
        }

        if (hmax > epsilon){    // Exclusively used for the LT_epsilon variant
            delete[] ind_1;
            delete[] ind_2;
            delete[] prior;
            delete[] unprior;
            return false;
        }
    }

    // Iterate over non-blacklisted, "normal" elements
    for (unsigned int i = 0; i < outer_loop_size; i++) {

        float hmin = INFINITY;

        for (unsigned int j = 0; j < m; j++) {

            float d;

            if(cache){
                d = flipped ? (*cache).cacheDistance(ind_2[j],ind_1[i]) : (*cache).cacheDistance(ind_1[i],ind_2[j]);
            }
            else{
                d = T::dist(t1.at(ind_1[i]), t2.at(ind_2[j]));
            }

            if (d < hmax) {
                hmin = INFINITY;   // when breaking the loop we have to reset hmin to preserve the integrity of the "hmin > hmax" check
                break;
            }

            if (d < hmin) {
                hmin = d;
            }

        }

        if (hmin > hmax && !(std::isinf(hmin))) {
            hmax = hmin;
        }

        if (hmax > epsilon){    // Exclusively used for the LT_epsilon variant
            delete[] ind_1;
            delete[] ind_2;
            delete[] prior;
            delete[] unprior;
            return false;
        }

    }

    // Iterate over unprioritized if not empty
    for(int i = 0; i < hd_ind.get_unprioritized()->size(); i++){
        // NN Search
        float hmin = INFINITY;

        for (unsigned int j = 0; j < m; j++) {

            float d;

            if(cache){
                d = flipped ? (*cache).cacheDistance(ind_2[j],unprior[i]) : (*cache).cacheDistance(unprior[i],ind_2[j]);
            }
            else{
                d = T::dist(t1.at(unprior[i]), t2.at(ind_2[j]));
            }

            if (d < hmax) {
                hmin = INFINITY;   // when breaking the loop we have to reset hmin to preserve the integrity of the "hmin > hmax" check
                break;
            }

            if (d < hmin) {
                hmin = d;
            }

        }

        if (hmin > hmax && !(std::isinf(hmin))) {
            hmax = hmin;
        }

        if (hmax > epsilon){    // Exclusively used for the LT_epsilon variant
            delete[] ind_1;
            delete[] ind_2;
            delete[] prior;
            delete[] unprior;
            return false;
        }
    }

    delete[] ind_1;
    delete[] ind_2;
    delete[] prior;
    delete[] unprior;

    return (hmax <= epsilon);

}

template<class T>
bool modified_hd_LT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float epsilon, HausdorffIndexManager& hd_ind, std::shared_ptr<DistCache<T>> cache = NULL){

    if(hd_ind.nr_blacklisted_i != t1.length()){
        bool dist_1 = modified_dirHd_LT_epsilon(t1, t2, false, epsilon, hd_ind, cache);

        if(!dist_1){
            return false;
        }

    }

    if(hd_ind.nr_blacklisted_j != t2.length()){
        bool dist_2 = modified_dirHd_LT_epsilon(t2, t1, true, epsilon, hd_ind, cache);

        return dist_2;

    }

    return true;

}