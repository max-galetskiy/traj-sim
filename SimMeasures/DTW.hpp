#pragma once
#include <algorithm>
#include "../TrajModels/Trajectory.hpp"
#include "../DistCache/DistCache.hpp"

template<class T>
float dynamic_time_warping(Trajectory<T>& t1, Trajectory<T>& t2, std::shared_ptr<DistCache<T>> cache = nullptr) {

	const unsigned int n = t1.length();
	const unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return (m == n) ? 0 : INFINITY;
    }

	auto* dtw = new float[n * m];

	for (unsigned int i = 0; i < n; i++) {

		for (unsigned int j = 0; j < m; j++) {

			float d = cache ? (*cache).cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));

			if (i == 0 && j == 0) {
				dtw[i * m + j] = d;																// dtw[i*m + j] === dtw[i][j]
			}
			else if (i == 0) {
				dtw[i * m + j] = d + dtw[i * m + (j - 1)];
			}
			else if (j == 0) {
				dtw[i * m + j] = d + dtw[(i - 1) * m + j];
			}
			else {
				dtw[i * m + j] = d + std::min({ dtw[i * m + (j - 1)], dtw[(i - 1) * m + j], dtw[(i - 1) * m + (j - 1)] });
			}
		}

	}

	float dist = dtw[(n - 1) * m + (m - 1)];
	delete[] dtw;
	return dist;

}

template<class T>
bool dtw_LT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float epsilon, std::shared_ptr<DistCache<T>> cache = nullptr){

	const unsigned int n = t1.length();
	const unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return (m == n);
    }

	auto* dtw = new float[n * m];

	bool curr_diag = true;				// True if all values are > epsilon
	bool prev_diag = false;

	// Traverse Matrix in diagonal fashion 
	for (int k = 0; k < (n + m - 1); k++){

		int i;
		int j; 

		// Early Break
		if(prev_diag && curr_diag){
            delete[] dtw;
            return false;
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
				dtw[i * m + j] = d;																
			}
			else if (i == 0) {
				dtw[i * m + j] = d + dtw[i * m + (j - 1)];
			}
			else if (j == 0) {
				dtw[i * m + j] = d + dtw[(i - 1) * m + j];
			}
			else {
				dtw[i * m + j] = d + std::min({ dtw[i * m + (j - 1)], dtw[(i - 1) * m + j], dtw[(i - 1) * m + (j - 1)] });
			}

			if(dtw[i * m + j] < epsilon){
				curr_diag = false;
			}

			i--;
			j++;
		}

	}

	float dist = dtw[(n - 1) * m + (m - 1)];
	delete[] dtw;
	return (dist < epsilon);

}