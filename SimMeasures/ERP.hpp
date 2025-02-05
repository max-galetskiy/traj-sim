#pragma once
#include <algorithm>
#include "../TrajModels/Trajectory.hpp"
#include "../DistCache/DistCache.hpp"

template<class T>
float edit_distance_real_penalty(Trajectory<T>& t1, Trajectory<T>& t2, T ref_point, std::shared_ptr<DistCache<T>> cache = nullptr) {

	const unsigned int n = t1.length();
	const unsigned int m = t2.length();

	auto* erp_empt = new float[n + m];		// Stores values for cases where one trajectory is empty

	// Fill erp_empt
	for (unsigned int i = 0; i < n; i++){
		float d = T::dist(t1.at(i), ref_point);
		erp_empt[i] = (i == 0) ? d : d + erp_empt[i - 1];
 	}
	for (unsigned int j = 0; j < m; j++){					// Values for i == 0 are put in the same array as values for j == 0 to allocate memory in one block
		float d = T::dist(ref_point, t2.at(j));
		erp_empt[n + j] = (j == 0) ? d : d + erp_empt[n + j - 1];
	}

    if(m == 0 || n == 0){
        if( m == n){
            delete[]erp_empt;
            return 0;
        }else{
            float distance = (m == 0) ? erp_empt[n - 1] : erp_empt[n + m - 1];
            delete[] erp_empt;
            return distance;
        }
    }

    auto* erp = new float[n * m];

    // Fill erp
	for (unsigned int i = 0; i < n; i++) {

		for (unsigned int j = 0; j < m; j++) {

				float d = cache ? (*cache).cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));

				float left = (j == 0) ? erp_empt[i] : erp[i * m + (j - 1)];

				float bottom = (i == 0) ? erp_empt[n + j] : erp[(i - 1) * m + j];

				float diag; 

				if((i > 0) && (j > 0)){
					diag = erp[(i - 1) * m + (j - 1)];
				}
				else if ((i == 0) && (j == 0)){
					diag = 0;
				}
				else{
					diag = (i == 0) ? erp_empt[n + (j - 1)] : erp_empt[i - 1];
				}

				erp[i * m + j] = std::min({ left + T::dist(ref_point, t2.at(j)), bottom + T::dist(t1.at(i), ref_point), diag + d});

		}
	}

	float ed = erp[(n - 1) * m + (m - 1)];
	delete[] erp;
	delete[] erp_empt;
	return ed;

}

template<class T>
bool erp_LT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, T ref_point, float epsilon, std::shared_ptr<DistCache<T>> cache = nullptr){

	const unsigned int n = t1.length();
	const unsigned int m = t2.length();

	auto* erp_empt = new float[n + m];		// Stores values for cases where one trajectory is empty

	// Fill erp_empt
	for (unsigned int i = 0; i < n; i++){
		float d = T::dist(t1.at(i), ref_point);
		erp_empt[i] = (i == 0) ? d : d + erp_empt[i - 1];
 	}
	for (unsigned int j = 0; j < m; j++){					// Values for i == 0 are put in the same array as values for j == 0 to allocate memory in one block
		float d = T::dist(ref_point, t2.at(j));
		erp_empt[n + j] = (j == 0) ? d : d + erp_empt[n + j - 1];
	}

    auto* erp = new float[n * m];

    bool curr_diag = true;				// True if all values are > epsilon
	bool prev_diag = false;

	// Traverse Matrix in diagonal fashion 
	for (int k = 0; k < (n + m - 1); k++){

		int i;
		int j; 

		// Early Break
		if(prev_diag && curr_diag){
            delete[] erp_empt;
            delete[] erp;
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

			float d = cache ? (*cache).cacheDistance(i, j) : T::dist(t1.at(i), t2.at(j));

			float left = (j == 0) ? erp_empt[i] : erp[i * m + (j - 1)];

			float bottom = (i == 0) ? erp_empt[n + j] : erp[(i - 1) * m + j];

			float diag;

			if ((i > 0) && (j > 0)){
				diag = erp[(i - 1) * m + (j - 1)];
			}
			else if ((i == 0) && (j == 0)){
				diag = 0;
			}
			else{
				diag = (i == 0) ? erp_empt[n + (j - 1)] : erp_empt[i - 1];
			}

			erp[i * m + j] = std::min({left + T::dist(ref_point, t2.at(j)), bottom + T::dist(t1.at(i), ref_point), diag + d});

			if (erp[i * m + j] < epsilon){
				curr_diag = false;
			}

			i--;
			j++;
		}

	}

	float ed = erp[(n - 1) * m + (m - 1)];
	delete[] erp;
    delete[] erp_empt;
    return (ed < epsilon);
}