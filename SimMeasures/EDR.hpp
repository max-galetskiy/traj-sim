#pragma once
#include <algorithm>
#include "../TrajModels/Trajectory.hpp"
#include "../DistCache/DistCache.hpp"

template<class T>
float edit_distance_real_sequence(Trajectory<T>& t1, Trajectory<T>& t2, float tau, std::shared_ptr<DistCache<T>> cache = nullptr) {

	const unsigned int n = t1.length();
	const unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return std::max(m,n); // returns 0 if both trajectories are empty, otherwise returns the length of the other trajectory
    }

	auto* edr = new float[n * m];

	for (unsigned int i = 0; i < n; i++) {

		for (unsigned int j = 0; j < m; j++) {

			float d = cache ? (*cache).cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));
			float cost   = (d < tau) ? 0 : 1;

			if (i == 0 && j == 0) {
				edr[i * m + j] = cost;														// edr[i*m + j] === edr[i][j]
			}
			else {

				float left   = (j > 0) ? edr[i * m + (j - 1)] : (i + 1);  // when one trajectory is empty, we fill in the length of the other
				float bottom = (i > 0) ? edr[(i - 1) * m + j] : (j + 1);
				
				float diag;

				if((i > 0) && (j > 0)){
					diag = edr[(i - 1) * m + (j - 1)];
				}
				else{
					diag = (i == 0) ? j : i;
				}

				edr[i * m + j] = std::min({ left + 1, bottom + 1, diag + cost });

			}
		}

	}

	float ed = edr[(n - 1) * m + (m - 1)];
	delete[] edr;
	return ed/std::max(m,n);

}


template<class T>
bool edr_LT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float tau, float epsilon, std::shared_ptr<DistCache<T>> cache = nullptr){

	const unsigned int n = t1.length();
	const unsigned int m = t2.length();


	auto* edr = new float[n * m];

	bool curr_diag = true;				// True if all values are > epsilon
	bool prev_diag = false;

	// Traverse Matrix in diagonal fashion 
	for (int k = 0; k < (n + m - 1); k++){

		int i;
		int j; 

		// Early Break
		if(prev_diag && curr_diag){
            delete[] edr;
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
			float cost   = (d < tau) ? 0 : 1;

			if (i == 0 && j == 0) {
				edr[i * m + j] = cost;														
			}
			else {

				float left   = (j > 0) ? edr[i * m + (j - 1)] : (i + 1);  // when one trajectory is empty, we fill in the length of the other
				float bottom = (i > 0) ? edr[(i - 1) * m + j] : (j + 1);
				
				float diag;

				if((i > 0) && (j > 0)){
					diag = edr[(i - 1) * m + (j - 1)];
				}
				else{
					diag = (i == 0) ? j : i;
				}

				edr[i * m + j] = std::min({ left + 1, bottom + 1, diag + cost });

			}

			if(edr[i * m + j] < epsilon){
				curr_diag = false;
			}


			i--;
			j++;
		}

	}

	float ed = edr[(n - 1) * m + (m - 1)];
	delete[] edr;
	return (ed < epsilon);

}