#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include "../TrajModels/Trajectory.hpp"
#include "../DistCache/DistCache.hpp"

template<class T>
float longest_common_subsequence(Trajectory<T>& t1, Trajectory<T>& t2, float tau, std::shared_ptr<DistCache<T>> cache = nullptr) {

	const unsigned int n = t1.length();
	const unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return 0;
    }

	auto* lcss = new float[n * m];

	for (unsigned int i = 0; i < n; i++) {

		for (unsigned int j = 0; j < m; j++) {

			float d = cache ? (*cache).cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));

			if (i == 0 && j == 0) {
				lcss[i * m + j] = (d < tau) ? 1 : 0;														// lcss[i*m + j] === lcss[i][j]
			}
			else if (i == 0) {
				lcss[i * m + j] = (d < tau) ? 1 : lcss[i * m + (j - 1)];
			}
			else if (j == 0) {
				lcss[i * m + j] = (d < tau) ? 1 : lcss[(i - 1) * m + j];
			}
			else {

				if (d < tau) {
					lcss[i * m + j] = 1 + lcss[(i - 1) * m + (j - 1)];
				}
				else {
					lcss[i * m + j] = std::max(lcss[(i - 1) * m + j], lcss[i * m + (j - 1)]);
				}
			}
		}

	}

	float len = lcss[(n - 1) * m + (m - 1)];
	delete[] lcss;
	return len/std::min(m,n);

}

/*
	Computes the restricted LCSS. Should not be called directly
*/
template<class T>
float internal_restricted_lcss(Trajectory<T>& t1, Trajectory<T>& t2, float tau, int delta, std::shared_ptr<DistCache<T>> cache = nullptr, float epsilon = -1){

    const unsigned int n = t1.length();
    const unsigned int m = t2.length();

    if(m == 0 || n == 0){
        return 0;
    }

    auto* lcss = new float[n * m];
    float max_value = -1;               // We store the maximum value so that we don't have to bother with tracking which columns we used

    for(unsigned int i = 0; i < n; i++){


        unsigned int min_col = std::max<int>(0, i - delta);
        unsigned int max_col = std::min<int>(m - 1, i + delta);

        bool row_greater_epsilon = true; // Used for the decision problem variant

        for(unsigned int j = min_col; j <= max_col; j++){

            float d = cache ? (*cache).cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));

            if (i == 0 && j == 0) {
                lcss[i * m + j] = (d < tau) ? 1 : 0;
            }
            else if (i == 0 && j > min_col) {
                lcss[i * m + j] = (d < tau) ? 1 : lcss[i * m + (j - 1)];
            }
            else if (j == 0) {
                lcss[i * m + j] = (d < tau) ? 1 : lcss[(i - 1) * m + j];
            }
            else {

                if (d < tau) {
                    lcss[i * m + j] = 1 + lcss[(i - 1) * m + (j - 1)];
                }
                else {

                    // handles cases where no left/top neighbour exists
                    float left = (j > min_col) ? lcss[i * m + (j - 1)] : 0;
                    float top  = (j < max_col || min_col == max_col) ? lcss[(i - 1) * m + j] : 0;               // Additional clause for the case that n > m + delta

                    lcss[i * m + j] = std::max(left,top);
                }
            }

            if(lcss[i * m + j] > max_value){
                max_value = lcss[i * m + j];
            }

            if(epsilon >= 0 && lcss[i * m + j] <= epsilon){
                row_greater_epsilon = false;
            }

        }

        if(epsilon >= 0 && row_greater_epsilon){
            delete[] lcss;
            return -1; // Value can be identified as early break
        }
    }

    delete[] lcss;
    return max_value;
}

template<class T>
float restricted_lcss(Trajectory<T>& t1, Trajectory<T>& t2, float tau, int delta, std::shared_ptr<DistCache<T>> cache = NULL){
    return internal_restricted_lcss(t1,t2,tau, delta, cache);
}

template<class T>
bool lcss_GT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float tau, float epsilon, std::shared_ptr<DistCache<T>> cache = NULL){

	const unsigned int n = t1.length();
	const unsigned int m = t2.length();

    if(n == 0 ||m == 0){
        return (epsilon == 0);
    }


    float* lcss = new float[n * m];

	bool curr_diag = true;				// True if all values are > epsilon
	bool prev_diag = false;

	// Traverse Matrix in diagonal fashion 
	for (int k = 0; k < (n + m - 1); k++){

		int i;
		int j; 

		// Early Break
		if(prev_diag && curr_diag){
            delete[] lcss;
			return true;
		}

		if(k < n){   // for diagonals before (and including) main diagonal
			i = k;
			j = 0;
		}
		else{
			i = n - 1;
			j = k - (n - 1);
		}

		prev_diag = curr_diag;
		curr_diag = true;

		// Traverse diagonal
		while(i >= 0 && j < m){
			
			float d = cache ? (*cache).cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));

			if (i == 0 && j == 0) {
				lcss[i * m + j] = (d < tau) ? 1 : 0;														
			}
			else if (i == 0) {
				lcss[i * m + j] = (d < tau) ? 1 : lcss[i * m + (j - 1)];
			}
			else if (j == 0) {
				lcss[i * m + j] = (d < tau) ? 1 : lcss[(i - 1) * m + j];
			}
			else {

				if (d < tau) {
					lcss[i * m + j] = 1 + lcss[(i - 1) * m + (j - 1)];
				}
				else {
					lcss[i * m + j] = std::max(lcss[(i - 1) * m + j], lcss[i * m + (j - 1)]);
				}
			}

			
			if(lcss[i * m + j] < epsilon){
				curr_diag = false;
			}

			i--;
			j++;
		}

	}

	float len = lcss[(n - 1) * m + (m - 1)];
	delete[] lcss;
	return (len >= epsilon);
}

template<class T>
bool restricted_lcss_GT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float tau, int delta, float epsilon = -1, std::shared_ptr<DistCache<T>> cache = NULL){

    float dist = internal_restricted_lcss(t1,t2,tau,delta,cache,epsilon);

    if(dist < 0){ // Signifies early break
        return true;
    }

    return (dist >= epsilon);
}