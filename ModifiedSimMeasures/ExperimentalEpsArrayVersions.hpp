#include <cmath>
#include <algorithm>
#include "../TrajModels/Trajectory.hpp"
#include "../DistCache/DistCache.hpp"

template<class T>
bool dFr_LT_epsilon_epsarray(Trajectory<T>& t1, Trajectory<T>& t2, float fr_eps, float hd_eps, float tau, std::shared_ptr<DistCache<T>> cache = NULL) {

    const unsigned int n = t1.length();
    const unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return (m == n);
    }

    float *dFr = new float[n*m];
    bool * epsarr_hd = new bool[n * m];
    bool * epsarr_tau = new bool[n * m];

    bool curr_diag = true;				// True if all values are > epsilon
    bool prev_diag = false;

    // Traverse Matrix in diagonal fashion
    for (int k = 0; k < (n + m - 1); k++){

        int i;
        int j;

        // Early Break
        if(prev_diag && curr_diag){
            delete[] dFr;
            delete[] epsarr_hd;
            delete[] epsarr_tau;
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

            epsarr_hd[i * m + j] = d <= hd_eps;
            epsarr_tau[i * m + j] = d <= tau;

            if (i == 0 && j == 0) {
                dFr[i * m + j] = d;
            }
            else if (i == 0) {
                dFr[i * m + j] = std::max(d, dFr[i * m + (j - 1)]);
            }
            else if (j == 0) {
                dFr[i * m + j] = std::max(d, dFr[(i - 1) * m + j]);
            }
            else {
                dFr[i * m + j] = std::max(d, std::min({ dFr[i * m + (j - 1)], dFr[(i - 1) * m + j], dFr[(i - 1) * m + (j - 1)] }));
            }

            if(dFr[i * m + j] < fr_eps){
                curr_diag = false;
            }

            i--;
            j++;
        }

    }

    float disc_fr = dFr[(n - 1) * m + (m - 1)];
    delete[] dFr;
    delete[] epsarr_hd;
    delete[] epsarr_tau;
    return (disc_fr < fr_eps);

}

template<class T>
bool lcss_GT_epsilon_epsarray(Trajectory<T>& t1, Trajectory<T>& t2, float fr_eps, float hd_eps, float tau, float lcss_eps, std::shared_ptr<DistCache<T>> cache = NULL){

    const unsigned int n = t1.length();
    const unsigned int m = t2.length();

    if(n == 0 ||m == 0){
        return (lcss_eps == 0);
    }

    float* lcss = new float[n * m];
    bool * epsarr_hd = new bool[n * m];
    bool * epsarr_fr = new bool[n * m];


    bool curr_diag = true;				// True if all values are > epsilon
    bool prev_diag = false;

    // Traverse Matrix in diagonal fashion
    for (int k = 0; k < (n + m - 1); k++){

        int i;
        int j;

        // Early Break
        if(prev_diag && curr_diag){
            delete[] lcss;
            delete[] epsarr_fr;
            delete[] epsarr_hd;
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

            epsarr_fr[i * m + j] = d <= fr_eps;
            epsarr_hd[i * m + j] = d <= hd_eps;

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


            if(lcss[i * m + j] < lcss_eps){
                curr_diag = false;
            }

            i--;
            j++;
        }

    }

    float len = lcss[(n - 1) * m + (m - 1)];
    delete[] lcss;
    delete[] epsarr_fr;
    delete[] epsarr_hd;
    return (len >= lcss_eps);

}

// template<class T>
// bool hd_LT_epsilon_epsarray(Trajectory<T>& t1, Trajectory<T>& t2, float fr_eps, float hd_eps, float tau, std::shared_ptr<DistCache<T>> cache = NULL)