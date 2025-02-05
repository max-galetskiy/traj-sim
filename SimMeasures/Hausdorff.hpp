#pragma once
#include <cmath>
#include <algorithm>
#include "../TrajModels/Trajectory.hpp"
#include "../DistCache/DistCache.hpp"
#include "../Utils/Utilities.hpp"

/*
    Computes the directed Hausdorff distance. This function should not be called directly!
*/
template<class T>
float dirHd(Trajectory<T>& t1, Trajectory<T>& t2, bool flipped, std::shared_ptr<DistCache<T>> cache = nullptr, float epsilon = -1) {

	// Missing: Remove Intersection

	const unsigned int n = t1.length();
	const unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return (m == n) ? 0 : INFINITY;
    }

	// Reorder Points
	int* ind_1 = new int[n];
	int* ind_2 = new int[m];
	inside_out_shuffle(ind_1,n);
	inside_out_shuffle(ind_2,m);

	// Calculate directed Hausdorff

	float hmax = 0;

	for (unsigned int i = 0; i < n; i++) {

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

        if (epsilon >= 0 && hmax > epsilon){    // Exclusively used for the LT_epsilon variant
            delete[] ind_1;
            delete[] ind_2;
            return -1;
        }

	}

	delete[] ind_1;
	delete[] ind_2;

	return hmax;

}

/*
    Computes the Hausdorff distance following the algorithm described by Taha and Hanbury
 */
template<class T>
float hausdorff(Trajectory<T>& t1, Trajectory<T>& t2, std::shared_ptr<DistCache<T>> cache = nullptr) {

    return std::max(dirHd(t1, t2, false, cache), dirHd(t2, t1, true, cache));

}

template<class T>
bool hd_LT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float epsilon, std::shared_ptr<DistCache<T>> cache = nullptr){

    float dist_1 = dirHd(t1, t2, false, cache, epsilon);

    if(dist_1 == -1){
        return false;
    }

    float dist_2 = dirHd(t2, t1, true, cache, epsilon);

    if(dist_2 == -1){
        return false;
    }

    return true;

}