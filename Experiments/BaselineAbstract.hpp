#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <array>
#include <utility>

#include "AbstractMeasure.hpp"

#include "../SimMeasures/Hausdorff.hpp"
#include "../SimMeasures/Frechet.hpp"
#include "../SimMeasures/DiscreteFrechet.hpp"
#include "../SimMeasures/DTW.hpp"
#include "../SimMeasures/LCSS.hpp"
#include "../SimMeasures/EDR.hpp"
#include "../SimMeasures/ERP.hpp"

#include "../ModifiedSimMeasures/ExperimentalEpsArrayVersions.hpp"

enum MEASURE {Frechet = 0, dFr, DTW, EDR, ERP, Hd, LCSS, rLCSS};

/*
    Note: Since I have removed the associated dataset, this code will not run (but a suitable call to load a dataset can be inserted quite easily if someone would want that)
    Given a trajectory model and a data loader, this computes and logs a couple of baselines
 */
template<class T, class DataLoader>
class BaselineAbstract {

    protected:

        inline static const int LW = 20; // whitespace width for log
        inline static const int CW = 35; // whitespace width for cout

        inline static std::shared_ptr<std::shared_ptr<Trajectory<T>>[]> trajectories;
        inline static std::fstream log;
        inline static std::string log_path = "../Experiments/baseline_logs/";

        inline static const std::string measure_string[] = {"Frechet", "Discrete Frechet", "Dynamic Time Warping", "Edit Distance on Real Sequences",
                                                        "Edit Distance with Real Penalty", "Hausdorff", "Longest Common Subsequence",
                                                        "Restricted Longest Common Subsequence"};

        // Log Functions
        static void init_log(std::string filename, std::fstream& temp_log);
        static void write_line(float value, std::fstream& temp_log);
        static void save_log(std::string filename, std::fstream& temp_log);

        // Formatting Function
        static std::string format_measure(MEASURE measure, T ref_point, float dist_threshold = -1,int ind_threshold = 0);

        // Computation Functions
        static double compute_measure(MEASURE measure, unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point);
        static double decide_measure(MEASURE measure, unsigned int trajectory_cardinality, float dist_threshold,int ind_threshold,T ref_point, float epsilon);

    public:
        static void set_log_path(std::string log_path);

        // Full Experiments
        static void compute_all_separately(unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point, bool write_file, bool ignore_frechet = false);
        static void decide_all_separately(unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point, std::array<float,8> epsilons, bool write_file);
        static void compute_all_separately_full_cache(unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point, bool write_file, bool ignore_frechet = false);
        static void decide_all_separately_full_cache(unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point, std::array<float,8> epsilons, bool write_file);

        static void compare_eps_array_versions(unsigned int trajectory_cardinality, float fr_eps, float hd_eps, float lcss_eps, float tau);
};


/*
    ++++++++++++++++++++++++++++++++
       COMPLETE BASELINES
    ++++++++++++++++++++++++++++++++
 */

template<class T, class DataLoader>
void BaselineAbstract<T, DataLoader>::compute_all_separately(unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point, bool write_file, bool ignore_frechet) {

    if(write_file){
        log = std::fstream(log_path + "baseline_1.txt", std::fstream::out | std::fstream::in | std::fstream::trunc);
        if(!log){
            std::throw_with_nested(std::runtime_error("Could not create baseline_1.txt!"));
            return;
        }

        if (!ignore_frechet) {log << TAB << "Fr";}
        log << TAB << "dFr" << TAB << "DTW" << TAB << "Hd" << TAB << "EDR" << TAB << "ERP" << TAB << "LCSS" << TAB << "resLCSS" << std::endl;
    }

    std::cout << "Starting Baseline 1: Compute all Measures separately" << std::endl;

    std::cout << "Loading Dataset...";
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Running Tests..." << std::endl;

    double frechet_runtime;
    if (!ignore_frechet) {frechet_runtime = compute_measure(Frechet, n, dist_threshold, ind_threshold,ref_point);}

    double dfr_runtime = compute_measure(dFr, n, dist_threshold, ind_threshold, ref_point);
    double dtw_runtime = compute_measure(DTW, n, dist_threshold, ind_threshold, ref_point);
    double hd_runtime = compute_measure(Hd, n, dist_threshold, ind_threshold, ref_point);
    double edr_runtime = compute_measure(EDR, n, dist_threshold, ind_threshold, ref_point);
    double erp_runtime = compute_measure(ERP, n, dist_threshold, ind_threshold, ref_point);
    double lcss_runtime = compute_measure(LCSS, n, dist_threshold, ind_threshold, ref_point);
    double rlcss_runtime = compute_measure(rLCSS, n, dist_threshold, ind_threshold, ref_point);

    std::cout << "Tests Complete!" << std::endl;

    std::cout << "Runtimes (microseconds): " << std::endl;
    std::cout << COUT_TAB << "Measure" << COUT_TAB << "Total" << COUT_TAB << "Average" << std::endl;
    if (!ignore_frechet) {std::cout << COUT_TAB << "Frechet" << COUT_TAB << frechet_runtime << COUT_TAB << frechet_runtime/(n* (n+1) * 0.5) << std::endl;}
    std::cout << COUT_TAB << "Discrete Frechet" << COUT_TAB << dfr_runtime << COUT_TAB << dfr_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Dynamic Time Warping" << COUT_TAB << dtw_runtime << COUT_TAB << dtw_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Hausdorff" << COUT_TAB << hd_runtime << COUT_TAB << hd_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Edit Distance on Real Sequences" << COUT_TAB << edr_runtime << COUT_TAB << edr_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Edit Distance with Real Penalty" << COUT_TAB << erp_runtime << COUT_TAB << erp_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Longest Common Subsequence" << COUT_TAB << lcss_runtime << COUT_TAB << lcss_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Restricted LCSS" << COUT_TAB << rlcss_runtime << COUT_TAB << rlcss_runtime/(n* (n+1) * 0.5) << std::endl;

    log.close();
}

template<class T, class DataLoader>
void BaselineAbstract<T, DataLoader>::decide_all_separately(unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point, std::array<float, 8> epsilons, bool write_file) {

    if(write_file){
        log = std::fstream(log_path + "baseline_3.txt", std::fstream::out | std::fstream::in | std::fstream::trunc);
        if(!log){
            std::throw_with_nested(std::runtime_error("Could not create baseline_3.txt!"));
            return;
        }

        log << TAB << "Fr <= " << epsilons[0] << TAB << "dFr <= " << epsilons[1] << TAB << "DTW <= " << epsilons[2] << TAB << "Hd <= " << epsilons[3]
            << TAB <<"EDR <= " << epsilons[4] << TAB << "ERP <= " << epsilons[5] << TAB << "LCSS >= " << epsilons[6] << TAB << "resLCSS >= " << epsilons[7] << std::endl;
    }

    std::cout << "Starting Baseline 3: Decide all Measures separately" << std::endl;

    std::cout << "Loading Dataset...";
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Running Tests..." << std::endl;

    double frechet_runtime = decide_measure(Frechet, n, dist_threshold, ind_threshold, ref_point, epsilons[0]);
    double dfr_runtime = decide_measure(dFr, n, dist_threshold, ind_threshold, ref_point, epsilons[1]);
    double dtw_runtime = decide_measure(DTW, n, dist_threshold, ind_threshold, ref_point, epsilons[2]);
    double hd_runtime = decide_measure(Hd, n, dist_threshold, ind_threshold, ref_point, epsilons[3]);
    double edr_runtime = decide_measure(EDR, n, dist_threshold, ind_threshold, ref_point, epsilons[4]);
    double erp_runtime = decide_measure(ERP, n, dist_threshold, ind_threshold, ref_point,epsilons[5]);
    double lcss_runtime = decide_measure(LCSS, n, dist_threshold, ind_threshold, ref_point, epsilons[6]);
    double rlcss_runtime = decide_measure(rLCSS, n, dist_threshold, ind_threshold, ref_point, epsilons[7]);

    std::cout << "Tests Complete!" << std::endl;

    std::cout << "Runtimes (microseconds): " << std::endl;
    std::cout << COUT_TAB << "Measure" << COUT_TAB << "Total" << COUT_TAB << "Average" << std::endl;
    std::cout << COUT_TAB << "Frechet" << COUT_TAB << frechet_runtime << COUT_TAB << frechet_runtime/(n* (n+1) * 0.5)  << std::endl;
    std::cout << COUT_TAB << "Discrete Frechet" << COUT_TAB << dfr_runtime << COUT_TAB << dfr_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Dynamic Time Warping" << COUT_TAB << dtw_runtime << COUT_TAB << dtw_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Hausdorff" << COUT_TAB << hd_runtime << COUT_TAB << hd_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Edit Distance on Real Sequences" << COUT_TAB << edr_runtime << COUT_TAB << edr_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Edit Distance with Real Penalty" << COUT_TAB << erp_runtime << COUT_TAB << erp_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Longest Common Subsequence" << COUT_TAB << lcss_runtime << COUT_TAB << lcss_runtime/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Restricted LCSS" << COUT_TAB << rlcss_runtime << COUT_TAB << rlcss_runtime/(n* (n+1) * 0.5) << std::endl;

    log.close();
}

template<class T, class DataLoader>
void BaselineAbstract<T, DataLoader>::compute_all_separately_full_cache(unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point, bool write_file, bool ignore_frechet){

    // Load Dataset
    if(write_file){
        log = std::fstream(log_path + "baseline_2.txt", std::fstream::out | std::fstream::in | std::fstream::trunc);
        if(!log){
            std::throw_with_nested(std::runtime_error("Could not create baseline_2.txt!"));
            return;
        }

        if(!ignore_frechet) {log << TAB << "Fr";}
        log << TAB << "dFr" << TAB << "DTW" << TAB << "Hd" << TAB << "EDR" << TAB << "ERP" << TAB << "LCSS" << TAB << "resLCSS" << std::endl;
    }

    std::cout << "Starting Baseline 2: Cache all distances and then compute measures separately" << std::endl;

    std::cout << "Loading Dataset...";
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Running Tests..." << std::endl;

    // Prepare Variables

    // Names
    std::string names[8];
    for(int i = 0; i < 8; i++){
        names[i] = format_measure(MEASURE(i), ref_point, dist_threshold, ind_threshold);
    }

    // Clocks
    double measure_runtimes[9];
    std::fill(measure_runtimes, measure_runtimes+9, 0);

    // Functions
    typedef Trajectory<T>& traj_ref;
    typedef std::shared_ptr<DistCache<T>> cache;
    float (*norm_measure_pntrs[4])(traj_ref , traj_ref, cache) = {frechet, discrete_frechet, dynamic_time_warping, hausdorff};
    // For edr, erp, lcss and rlcss no extra pointers are needed

    // Start Computation
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + trajectory_cardinality; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + trajectory_cardinality; j_it++) {

            // Compute Cache
            auto cache_start = std::chrono::high_resolution_clock::now();
            auto c = std::make_shared<DistCache<T>>(*i_it, *j_it, true);
            auto cache_end = std::chrono::high_resolution_clock::now();
            measure_runtimes[0] += std::chrono::duration_cast<std::chrono::microseconds>(cache_end - cache_start).count();

            // Compute all Measures using the cache

            // Measures without extra parameters
            for(int i = (!ignore_frechet) ? 0 : 1; i < 4; i++){
                auto start = std::chrono::high_resolution_clock::now();
                float d = norm_measure_pntrs[i](**i_it, **j_it, c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[i + 1] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << d;
                }
            }

            // EDR
            {
                auto start = std::chrono::high_resolution_clock::now();
                float d = edit_distance_real_sequence(**i_it, **j_it, dist_threshold, c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[5] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << d;
                }
            }

            // ERP
            {
                auto start = std::chrono::high_resolution_clock::now();
                float d = edit_distance_real_penalty(**i_it, **j_it, ref_point, dist_threshold, c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[6] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << d;
                }
            }

            // LCSS
            {
                auto start = std::chrono::high_resolution_clock::now();
                float d = longest_common_subsequence(**i_it, **j_it, dist_threshold, c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[7] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << d;
                }
            }

            // rLCSS
            {
                auto start = std::chrono::high_resolution_clock::now();
                float d = restricted_lcss(**i_it, **j_it, dist_threshold, ind_threshold, c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[8] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << d << std::endl;
                }
            }


        }
    }

    std::cout << "Tests Complete!" << std::endl;

    std::cout << "Runtimes (microseconds): " << std::endl;
    std::cout << COUT_TAB << "Measure" << COUT_TAB << "Total" << COUT_TAB << "Average" << std::endl;
    std::cout << COUT_TAB << "Cross-Product" << COUT_TAB << measure_runtimes[0] << COUT_TAB << measure_runtimes[0]/(n* (n+1) * 0.5)<< std::endl;
    if (!ignore_frechet) {std::cout << COUT_TAB << "Frechet" << COUT_TAB << measure_runtimes[1] << COUT_TAB << measure_runtimes[1]/(n* (n+1) * 0.5) << std::endl;}
    std::cout << COUT_TAB << "Discrete Frechet" << COUT_TAB << measure_runtimes[2] << COUT_TAB << measure_runtimes[2]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Dynamic Time Warping" << COUT_TAB << measure_runtimes[3] << COUT_TAB << measure_runtimes[3]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Hausdorff" << COUT_TAB << measure_runtimes[4] << COUT_TAB << measure_runtimes[4]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Edit Distance on Real Sequences" << COUT_TAB << measure_runtimes[5] << COUT_TAB << measure_runtimes[5]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Edit Distance with Real Penalty" << COUT_TAB << measure_runtimes[6] << COUT_TAB << measure_runtimes[6]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Longest Common Subsequence" << COUT_TAB << measure_runtimes[7] << COUT_TAB << measure_runtimes[7]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Restricted LCSS" << COUT_TAB << measure_runtimes[8] << COUT_TAB << measure_runtimes[8]/(n* (n+1) * 0.5) << std::endl;
}

template<class T, class DataLoader>
void BaselineAbstract<T, DataLoader>::decide_all_separately_full_cache(unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point, std::array<float,8> epsilons, bool write_file){

    // Load Dataset
    if(write_file){
        log = std::fstream(log_path + "baseline_4.txt", std::fstream::out | std::fstream::in | std::fstream::trunc);
        if(!log){
            std::throw_with_nested(std::runtime_error("Could not create baseline_4.txt!"));
            return;
        }

        log << TAB << "Fr" << TAB << "dFr" << TAB << "DTW" << TAB << "Hd" << TAB << "EDR" << TAB << "ERP" << TAB << "LCSS" << TAB << "resLCSS" << std::endl;
    }

    std::cout << "Starting Baseline 4: Cache all distances and then decide measures separately" << std::endl;

    std::cout << "Loading Dataset...";
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Running Tests..." << std::endl;

    // Prepare Variables

    // Names
    std::string names[8];
    for(int i = 0; i < 8; i++){
        names[i] = format_measure(MEASURE(i), ref_point, dist_threshold, ind_threshold);
    }

    // Clocks
    double measure_runtimes[9];
    std::fill(measure_runtimes, measure_runtimes+9, 0);

    // Functions
    typedef Trajectory<T>& traj_ref;
    typedef std::shared_ptr<DistCache<T>> cache;
    bool (*norm_measure_pntrs[4])(traj_ref , traj_ref, float, cache) = {frechet_LT_Epsilon, dFr_LT_epsilon,dtw_LT_epsilon,hd_LT_epsilon};
    // For edr, erp, lcss and rlcss no extra pointers are needed

    // Start Computation
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + trajectory_cardinality; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + trajectory_cardinality; j_it++) {

            // Compute Cache
            auto cache_start = std::chrono::high_resolution_clock::now();
            auto c = std::make_shared<DistCache<T>>(*i_it, *j_it, true);
            auto cache_end = std::chrono::high_resolution_clock::now();
            measure_runtimes[0] += std::chrono::duration_cast<std::chrono::microseconds>(cache_end - cache_start).count();

            // Compute all Measures using the cache

            // Measures without extra parameters
            for(int i = 0; i < 4; i++){
                auto start = std::chrono::high_resolution_clock::now();
                bool p = norm_measure_pntrs[i](**i_it, **j_it, epsilons[i], c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[i + 1] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << p;
                }
            }

            // EDR
            {
                auto start = std::chrono::high_resolution_clock::now();
                bool p = edr_LT_epsilon(**i_it, **j_it, dist_threshold, epsilons[4], c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[5] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << p;
                }
            }

            // ERP
            {
                auto start = std::chrono::high_resolution_clock::now();
                bool p = erp_LT_epsilon(**i_it, **j_it, ref_point, dist_threshold, epsilons[5], c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[6] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << p;
                }
            }

            // LCSS
            {
                auto start = std::chrono::high_resolution_clock::now();
                bool p = lcss_GT_epsilon(**i_it, **j_it, dist_threshold, epsilons[6], c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[7] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << p;
                }
            }

            // rLCSS
            {
                auto start = std::chrono::high_resolution_clock::now();
                bool p = restricted_lcss_GT_epsilon(**i_it, **j_it, dist_threshold, ind_threshold, epsilons[7], c);
                auto end = std::chrono::high_resolution_clock::now();
                measure_runtimes[8] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                if(write_file){
                    log << TAB << p << std::endl;
                }
            }


        }
    }

    std::cout << "Tests Complete!" << std::endl;

    std::cout << "Runtimes (microseconds): " << std::endl;
    std::cout << COUT_TAB << "Measure" << COUT_TAB << "Total" << COUT_TAB << "Average" << std::endl;
    std::cout << COUT_TAB << "Cross-Product" << COUT_TAB << measure_runtimes[0] << COUT_TAB << measure_runtimes[0]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Frechet" << COUT_TAB << measure_runtimes[1] << COUT_TAB << measure_runtimes[1]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Discrete Frechet" << COUT_TAB << measure_runtimes[2] << COUT_TAB << measure_runtimes[2]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Dynamic Time Warping" << COUT_TAB << measure_runtimes[3] << COUT_TAB << measure_runtimes[3]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Hausdorff" << COUT_TAB << measure_runtimes[4] << COUT_TAB << measure_runtimes[4]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Edit Distance on Real Sequences" << COUT_TAB << measure_runtimes[5] << COUT_TAB << measure_runtimes[5]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Edit Distance with Real Penalty" << COUT_TAB << measure_runtimes[6] << COUT_TAB << measure_runtimes[6]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Longest Common Subsequence" << COUT_TAB << measure_runtimes[7] << COUT_TAB << measure_runtimes[7]/(n* (n+1) * 0.5) << std::endl;
    std::cout << COUT_TAB << "Restricted LCSS" << COUT_TAB << measure_runtimes[8] << COUT_TAB << measure_runtimes[8]/(n* (n+1) * 0.5) << std::endl;
}

template<class T, class DataLoader>
void BaselineAbstract<T, DataLoader>::compare_eps_array_versions(unsigned int trajectory_cardinality, float fr_eps,float hd_eps, float lcss_eps, float tau) {

    // Load Dataset
    std::cout << "Starting Experiment: Evaluate excess runtime from epsilon array computation" << std::endl;

    std::cout << "Loading Dataset...";
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Running Tests..." << std::endl;

    // Clocks
    double normal_runtimes[3] = {0, 0, 0};
    double epsarr_runtimes[3] = {0, 0, 0};
    std::chrono::system_clock::time_point start, end;

    int p = 0;

    for (auto i_it = trajectories.get(); i_it < trajectories.get() + trajectory_cardinality; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + trajectory_cardinality; j_it++) {

            p++;

            start = std::chrono::high_resolution_clock::now();
            dFr_LT_epsilon(**i_it,**j_it,fr_eps);
            end = std::chrono::high_resolution_clock::now();
            normal_runtimes[0] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            start = std::chrono::high_resolution_clock::now();
            lcss_GT_epsilon(**i_it,**j_it,tau,lcss_eps);
            end = std::chrono::high_resolution_clock::now();
            normal_runtimes[1] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            start = std::chrono::high_resolution_clock::now();
            dFr_LT_epsilon_epsarray(**i_it,**j_it,fr_eps, hd_eps, tau);
            end = std::chrono::high_resolution_clock::now();
            epsarr_runtimes[0] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            start = std::chrono::high_resolution_clock::now();
            lcss_GT_epsilon_epsarray(**i_it,**j_it,fr_eps,hd_eps,tau,lcss_eps);
            end = std::chrono::high_resolution_clock::now();
            epsarr_runtimes[1] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();


        }
    }

    std::cout << "Done!" << std::endl << std::endl << "Results: (Total RT in ms, average in microseconds)" << std::endl;
    std::cout << TAB << "Measure" << TAB << "RT" << TAB << "RT Avg." << TAB << "RT Eps." << TAB << "RT Eps Avg." << std::endl;
    std::cout << TAB << "Disc. Frechet" << TAB << normal_runtimes[0]/1000 << TAB << normal_runtimes[0]/p << TAB << epsarr_runtimes[0]/1000 << TAB << epsarr_runtimes[0]/p << std::endl;
    std::cout << TAB << "LCSS" << TAB << normal_runtimes[1]/1000 << TAB << normal_runtimes[1]/p << TAB << epsarr_runtimes[1]/1000 << TAB << epsarr_runtimes[1]/p << std::endl;

}
/*
    ++++++++++++++++++++++++++++++++
       SINGLE MEASURE COMPUTATION
    ++++++++++++++++++++++++++++++++
 */
template<class T, class DataLoader>
double BaselineAbstract<T, DataLoader>::compute_measure(MEASURE measure, unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point) {

    // Prepare runtime counter
    double total_runtime = 0;

    // Prepare log
    std::fstream temp_log;
    init_log("baseline_1", temp_log);

    // Prepare Console baseline_logs
    std::string name = format_measure(measure, ref_point, dist_threshold, ind_threshold);

    std::cout << "Computing " << name << "...";

    // Choose correct distance measure
    switch (measure) {
        case Frechet:
            AbstractMeasure<T>::setMeasure(frechet);
            break;
        case dFr:
            AbstractMeasure<T>::setMeasure(discrete_frechet);
            break;
        case DTW:
            AbstractMeasure<T>::setMeasure(dynamic_time_warping);
            break;
        case Hd:
            AbstractMeasure<T>::setMeasure(hausdorff);
            break;
        case EDR:
            AbstractMeasure<T>::setMeasure(edit_distance_real_sequence, dist_threshold);
            break;
        case ERP:
            AbstractMeasure<T>::setMeasure(edit_distance_real_penalty, dist_threshold, ref_point);
            break;
        case LCSS:
            AbstractMeasure<T>::setMeasure(longest_common_subsequence, dist_threshold);
            break;
        case rLCSS:
            AbstractMeasure<T>::setMeasure(restricted_lcss, dist_threshold, ind_threshold);
            break;
    }

    for (auto i_it = trajectories.get(); i_it < trajectories.get() + trajectory_cardinality; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + trajectory_cardinality; j_it++) {

            auto start = std::chrono::high_resolution_clock::now();
            float d = AbstractMeasure<T>::calc_dist(**i_it, **j_it, nullptr);
            auto end = std::chrono::high_resolution_clock::now();

            total_runtime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            if(temp_log){
                write_line(d,temp_log);
            }

        }
    }

    std::cout << "Complete!" << std::endl;

    // Overwrite changes to log
    if(temp_log){
        save_log("baseline_1", temp_log);
    }

    return total_runtime;
}

template<class T, class DataLoader>
double BaselineAbstract<T, DataLoader>::decide_measure(MEASURE measure, unsigned int trajectory_cardinality, float dist_threshold, int ind_threshold, T ref_point, float epsilon) {

    // Prepare runtime counter
    double total_runtime = 0;

    // Prepare log
    std::fstream temp_log;
    init_log("baseline_3", temp_log);

    // Prepare Console baseline_logs
    std::string name = format_measure(measure, ref_point, dist_threshold, ind_threshold);

    std::cout << "Deciding " << name << " [eps = " << epsilon << "]...";

    switch (measure) {
        case Frechet:
            AbstractMeasure<T>::setDecider(frechet_LT_Epsilon);
            break;
        case dFr:
            AbstractMeasure<T>::setDecider(dFr_LT_epsilon);
            break;
        case DTW:
            AbstractMeasure<T>::setDecider(dtw_LT_epsilon);
            break;
        case Hd:
            AbstractMeasure<T>::setDecider(hd_LT_epsilon);
            break;
        case EDR:
            AbstractMeasure<T>::setDecider(edr_LT_epsilon, dist_threshold);
            break;
        case ERP:
            AbstractMeasure<T>::setDecider(edr_LT_epsilon, dist_threshold);
            AbstractMeasure<T>::setDecider(erp_LT_epsilon, dist_threshold, ref_point);
            break;
        case LCSS:
            AbstractMeasure<T>::setDecider(lcss_GT_epsilon, dist_threshold);
            break;
        case rLCSS:
            AbstractMeasure<T>::setDecider(restricted_lcss_GT_epsilon, dist_threshold, ind_threshold);

    }

    for (auto i_it = trajectories.get(); i_it < trajectories.get() + trajectory_cardinality; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + trajectory_cardinality; j_it++) {
            auto start = std::chrono::high_resolution_clock::now();
            float d = AbstractMeasure<T>::decide_dist(**i_it, **j_it, epsilon, nullptr);
            auto end = std::chrono::high_resolution_clock::now();

            total_runtime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            if(temp_log){
                write_line(d,temp_log);
            }

        }
    }

    std::cout << "Complete!" << std::endl;

    // Overwrite changes to log
    if(temp_log){
        save_log("baseline_3", temp_log);
    }

    return total_runtime;
}

/*
    ++++++++++++++++++++++++++++++++
            FORMAT FUNCTION
    ++++++++++++++++++++++++++++++++
 */
template<class T, class DataLoader>
std::string BaselineAbstract<T, DataLoader>::format_measure(MEASURE measure, T ref_point, float dist_threshold, int ind_threshold) {

    std::string name = measure_string[measure];
    if(measure == EDR || measure == ERP || measure == LCSS || measure == rLCSS){
        name += " (tau = " + std::to_string(dist_threshold);

        if(measure == rLCSS){
            name += ", delta = " + std::to_string(ind_threshold) + ")";
        }
        if(measure == ERP){
            name += ", g = " + ref_point.to_string() + ")";
        }else{
            name += ")";
        }
    }
    return name;
}

/*
    ++++++++++++++++++++++++++++++++
            LOG FUNCTIONS
    ++++++++++++++++++++++++++++++++
 */
template<class T, class DataLoader>
void BaselineAbstract<T, DataLoader>::set_log_path(std::string log_path) {
    BaselineAbstract<T,DataLoader>::log_path = std::move(log_path);
}


template<class T, class DataLoader>
void BaselineAbstract<T, DataLoader>::init_log(std::string filename, std::fstream &temp_log) {
    // Prepare log
    if(log){

        log.clear();
        log.seekg(0);

        std::string path = log_path + filename + "_t.txt";
        temp_log = std::fstream(path, std::fstream::out); // create temporary log

        std::string line;               // skip first line
        std::getline(log,line);
        temp_log << line << std::endl;
    }
}


template<class T, class DataLoader>
void BaselineAbstract<T, DataLoader>::write_line(float value, std::fstream &temp_log) {

    std::string line;
    std::getline(log, line);

    if(line.empty()){
        temp_log << TAB << value << std::endl;
    }else{
        temp_log << line << TAB << value << std::endl;
    }
}


template<class T, class DataLoader>
void BaselineAbstract<T, DataLoader>::save_log(std::string filename, std::fstream &temp_log) {

    log.close();
    temp_log.close();
    std::string path = log_path + filename + ".txt";
    std::string temp_path = log_path + filename + "_t.txt";
    std::remove(path.c_str());
    std::rename(temp_path.c_str(), path.c_str());

    log.open(path, std::fstream::out | std::fstream ::in); // reopen log
}