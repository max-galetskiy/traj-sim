#pragma once
#include <memory>
#include <array>
#include <fstream>
#include <iostream>
#include <iomanip>


#include "../TrajModels/Trajectory.hpp"

#include "../SimMeasures/Frechet.hpp"
#include "../SimMeasures/Hausdorff.hpp"
#include "../SimMeasures/DTW.hpp"

#include "../ModifiedSimMeasures/CombinedSimilarityMeasures.hpp"
#include "../ModifiedSimMeasures/ExperimentalVersions.hpp"

/*
    Note: This entire file is mostly kept in here for documentation purposes. Since I have removed the associated dataset, this code will not run (but a suitable call to load a dataset can be inserted quite easily if someone would want that)
 */

template<class T, class DataLoader>
class CombinedMeasuresExperiments {

    private:
        inline static std::string log_path = "../Experiments/combined_measure_tests/";
        inline static const int LW = 20; // whitespace width for log
        inline static const int CW = 35; // whitespace width for cout

    public:

        static void compare_algorithms(unsigned int trajectory_cardinality, std::array<float,3> epsilons, bool write_file = false);

        static void pure_runtime(unsigned int trajectory_cardinality, std::array<float,3> epsilons);
        static void join_runtime(unsigned int trajectory_cardinality, std::array<float,3> epsilons, COMB_MEASURE combiner);
        static void decision_weighting(unsigned  int trajectory_cardinality, std::array<float, 3> epsilons, COMB_MEASURE   );
        static void fs_exploral(unsigned int trajectory_cardinality, std::array<float, 3> epsilons);
        static void selectivity(unsigned int trajectory_cardinality, std::array<float, 3> epsilons);

        static void pure_runtime_with_dfr(unsigned int trajectory_cardinality, std::array<float, 3> epsilons);
};

template<class T, class DataLoader>
void CombinedMeasuresExperiments<T, DataLoader>::compare_algorithms(unsigned int trajectory_cardinality, std::array<float, 3> epsilons, bool write_file) {

    std::cout << "Running Experiment: Comparison Between Independent Computation and Fr-Hd-DTW Algorithm" << std::endl;

    std::fstream comb_values_log;
    std::fstream comb_runtime_log;
    std::fstream base_values_log;
    std::fstream base_runtime_log;

    // Load Data
    std::cout << "Loading Dataset...";
    std::shared_ptr<std::shared_ptr<Trajectory<T>>[]> trajectories;
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Running Tests...";
    // Init Logs (if needed)
    if(write_file){
        base_values_log = std::fstream(log_path + "base_values.txt", std::fstream::out | std::fstream::in | std::fstream::trunc);
        comb_values_log = std::fstream(log_path + "comb_values.txt", std::fstream::out | std::fstream::in | std::fstream::trunc);
        base_runtime_log = std::fstream(log_path + "base_runtimes.txt", std::fstream::out | std::fstream::in | std::fstream::trunc);
        comb_runtime_log = std::fstream(log_path + "comb_runtimes.txt", std::fstream::out | std::fstream::in | std::fstream::trunc);

        if(!base_runtime_log || !base_values_log || !comb_values_log || !comb_runtime_log){
            std::throw_with_nested(std::runtime_error("Could not create neccessary log files!"));
            return;
        }

    }

    // Prepare Runtime Counters
    double base_runtimes[3] = {0,0,0};
    //double comb_runtimes[3] = {0,0,0};
    double comb_runtime = 0;

    // Iterate over trajectories
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {


            // Basic Computation
            bool base_fr = frechet_LT_Epsilon(**i_it, **j_it, epsilons[0]);
            bool base_hd = hd_LT_epsilon(**i_it, **j_it, epsilons[1]);
            bool base_dtw = dtw_LT_epsilon(**i_it,**j_it, epsilons[2]);


            // Using the Frechet-Hausdorff-DTW Algorithm
            auto comb_bools = frechet_hausdorff_dtw_algorithm(**i_it, **j_it, epsilons[0], epsilons[1], epsilons[2]);

            // Logging
            if(write_file){
                base_values_log << TAB << base_fr << TAB << base_hd << TAB << base_dtw << std::endl;
                comb_values_log << TAB << get<0>(comb_bools) << TAB << get<1>(comb_bools) << TAB << get<2>(comb_bools) << std::endl;
            }


        }
    }

    std::cout << "Done!" << std::endl;

    if(write_file){
        base_values_log.close();
        comb_values_log.close();
        base_runtime_log.close();
        comb_runtime_log.close();
    }

}


template<class T, class DataLoader>
void CombinedMeasuresExperiments<T, DataLoader>::pure_runtime(unsigned int trajectory_cardinality, std::array<float, 3> epsilons) {

    std::cout << "Running Experiment: Combined vs. Independent Computation - Pure Runtime" << std::endl;

    // Load Data
    std::cout << "Loading Dataset...";
    std::shared_ptr<std::shared_ptr<Trajectory<T>>[]> trajectories;
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Setup:" << std::endl;
    std::cout << "n = " << trajectory_cardinality << "; eps_fr = " << epsilons[0] << "; eps_dtw = " << epsilons[1] << "; eps_hd = " << epsilons[2] << std::endl;

    std::cout << "Running Tests...";

    // Setup Runtime Variables
    double base_runtime = 0, comb_runtime = 0, base_rt_best = INFINITY, comb_rt_best = INFINITY, base_rt_worst = 0, comb_rt_worst = 0;
    double base_fr = 0, base_hd = 0, base_dtw = 0, comb_fr = 0, comb_hd = 0, comb_dtw = 0;
    int compared_trajectories = 0;

    // Iterate over trajectories - Base Computation
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {

            double fr_rt = 0, hd_rt = 0, dtw_rt = 0;

            auto cache = std::make_shared<DistCache<T>>(*i_it, *j_it);

            auto start = std::chrono::high_resolution_clock::now();
            frechet_LT_Epsilon(**i_it, **j_it, epsilons[0], cache);
            auto end = std::chrono::high_resolution_clock::now();
            fr_rt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            base_fr += fr_rt;

            start = std::chrono::high_resolution_clock::now();
            hd_LT_epsilon(**i_it,**j_it,epsilons[2], cache);
            end = std::chrono::high_resolution_clock::now();
            hd_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            base_hd += hd_rt;

            start = std::chrono::high_resolution_clock::now();
            dtw_LT_epsilon(**i_it,**j_it,epsilons[1], cache);
            end = std::chrono::high_resolution_clock::now();
            dtw_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            base_dtw += dtw_rt;


            double rt = fr_rt + hd_rt/1000 + dtw_rt/1000;
            base_runtime += rt;

            if(rt < base_rt_best) {base_rt_best = rt;}
            if(rt > base_rt_worst) {base_rt_worst = rt;}

            compared_trajectories++;

        }
    }

    // Iterate over trajectories - Combined Computation
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {

            std::tuple<double, double, double> runtimes;

            runtimes = timed_fr_hd_dtw_algorithm(**i_it, **j_it, epsilons[0], epsilons[2], epsilons[1],true);
            comb_fr += get<0>(runtimes);
            comb_hd += get<1>(runtimes);
            comb_dtw += get<2>(runtimes);

            double rt = get<0>(runtimes) + get<1>(runtimes)/1000 + get<2>(runtimes)/1000;
            comb_runtime += rt;

            if(rt < comb_rt_best) {comb_rt_best = rt;}
            if(rt > comb_rt_worst) {comb_rt_worst = rt;}
        }
    }

    std::cout << "Done!" << std::endl;
    std::cout << "Trajectory pairs checked: " << compared_trajectories << std::endl;
    std::cout << std::endl;
    std::cout << "Results:" << std::endl << std::endl;

    std::cout << COUT_TAB << "RT. Indep." << COUT_TAB << "RT. Indep. Avg." << COUT_TAB << "RT. Indep. Best" << COUT_TAB << "RT. Indep. Worst" << COUT_TAB << "RT. Comb." << COUT_TAB << "RT. Comb. Avg." << COUT_TAB << "RT. Comb. Best" << COUT_TAB << "RT. Comb. Worst." << std::endl;
    std::cout << COUT_TAB << base_runtime << COUT_TAB << base_runtime/compared_trajectories << COUT_TAB << base_rt_best << COUT_TAB << base_rt_worst << COUT_TAB << comb_runtime << COUT_TAB << comb_runtime/compared_trajectories << COUT_TAB << comb_rt_best << COUT_TAB << comb_rt_worst << std::endl;
    std::cout << std::endl;
    std::cout << COUT_TAB << "Indep. Fr" << COUT_TAB << "Indep. Hd." << COUT_TAB << "Indep. DTW" << COUT_TAB << "Comb. Fr" << COUT_TAB << "Comb. Hd." << COUT_TAB << "Comb. DTW" << std::endl;
    std::cout << COUT_TAB << base_fr << COUT_TAB << base_hd/1000 << COUT_TAB << base_dtw/1000 << COUT_TAB << comb_fr << COUT_TAB << comb_hd/1000 << COUT_TAB << comb_dtw/1000 << std::endl;
    std::cout << std::endl << std::endl;

}

template<class T, class DataLoader>
void CombinedMeasuresExperiments<T, DataLoader>::decision_weighting(unsigned int trajectory_cardinality,std::array<float, 3> epsilons,COMB_MEASURE combiner) {

    std::cout << "Running Experiment: Impact of Positive vs. Negative Decision" << std::endl;

    // Load Data
    std::cout << "Loading Dataset...";
    std::shared_ptr<std::shared_ptr<Trajectory<T>>[]> trajectories;
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::string combiner_txt[3] = {"AND", "OR", "MAJORITY"};

    std::cout << "Setup:" << std::endl;
    std::cout << "n = " << trajectory_cardinality << "; eps_fr = " << epsilons[0] << "; eps_dtw = " << epsilons[1] << "; eps_hd = " << epsilons[2] << "; Combiner: " << combiner_txt[combiner] << std::endl;

    std::cout << "Running Tests...";


    // Setup Parameters
    double pos_rt = 0;
    double pos_best = INFINITY;
    double pos_worst = 0;
    double neg_rt = 0;
    double neg_best = INFINITY;
    double neg_worst = 0;
    int pos_trajectories = 0;
    int neg_trajectories = 0;

    // Iterate over trajectories
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {

            auto start = std::chrono::high_resolution_clock::now();
            bool d = fr_hd_dtw_similarity_join(**i_it, **j_it, epsilons[0], epsilons[2], epsilons[1], combiner, true);
            auto end = std::chrono::high_resolution_clock::now();

            double rt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            if(d){
                pos_trajectories++;
                pos_rt += rt;
                if(rt < pos_best){pos_best = rt;}
                if(rt > pos_worst){pos_worst = rt;}
            }else{
                neg_trajectories++;
                neg_rt += rt;
                if(rt < neg_best){neg_best = rt;}
                if(rt > neg_worst){neg_worst = rt;}
            }

        }
    }

    std::cout << "Done!" << std::endl;
    std::cout << "Positive Trajectory pairs: " << pos_trajectories << " || Negative Trajectory pairs: " << neg_trajectories << std::endl;
    std::cout << std::endl;
    std::cout << "Results:" << std::endl << std::endl;

    std::cout << COUT_TAB << "+ Best RT." << COUT_TAB << "+ Worst RT." << COUT_TAB << "+ Avg. RT." << COUT_TAB << "- Best RT." << COUT_TAB << "- Worst RT." << COUT_TAB << "- Avg. RT." << std::endl;
    std::cout << COUT_TAB << pos_best << COUT_TAB << pos_worst << COUT_TAB << pos_rt/pos_trajectories << COUT_TAB << neg_best << COUT_TAB << neg_worst << COUT_TAB << neg_rt/neg_trajectories << std::endl;
}

template<class T, class DataLoader>
void CombinedMeasuresExperiments<T, DataLoader>::join_runtime(unsigned int trajectory_cardinality, std::array<float, 3> epsilons, COMB_MEASURE combiner) {

    std::cout << "Running Experiment: Combined vs. Independent Computation - Similarity Join" << std::endl;

    // Load Data
    std::cout << "Loading Dataset...";
    std::shared_ptr<std::shared_ptr<Trajectory<T>>[]> trajectories;
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;
    std::string combiner_txt[3] = {"AND", "OR", "MAJORITY"};

    std::cout << "Setup:" << std::endl;
    std::cout << "n = " << trajectory_cardinality << "; eps_fr = " << epsilons[0] << "; eps_dtw = " << epsilons[1] << "; eps_hd = " << epsilons[2] << "; Combiner: " << combiner_txt[combiner] << std::endl;

    std::cout << "Running Tests...";

    // Setup Runtime Variables
    double base_best_case = 0, base_worst_case = 0, comb_rt = 0;
    double base_fr = 0, base_hd = 0, base_dtw = 0, comb_fr = 0, comb_hd = 0, comb_dtw = 0;
    int compared_trajectories = 0;

    // Iterate over trajectories - Base Computation
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {

            double fr_rt = 0, hd_rt = 0, dtw_rt = 0;
            bool fr, hd, dtw;
            double curr_best_case = 0, curr_worst_case = 0;

//            auto cache = std::make_shared<DistCache<T>>(*i_it, *j_it);

            auto start = std::chrono::high_resolution_clock::now();
            fr = frechet_LT_Epsilon(**i_it, **j_it, epsilons[0]);
            auto end = std::chrono::high_resolution_clock::now();
            fr_rt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            base_fr += fr_rt;

            start = std::chrono::high_resolution_clock::now();
            dtw = dtw_LT_epsilon(**i_it,**j_it,epsilons[1]);
            end = std::chrono::high_resolution_clock::now();
            dtw_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            base_dtw += dtw_rt;

            start = std::chrono::high_resolution_clock::now();
            hd = hd_LT_epsilon(**i_it,**j_it,epsilons[2]);
            end = std::chrono::high_resolution_clock::now();
            hd_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            base_hd += hd_rt;

            switch (combiner) {
                case AND:

                    if(fr && hd && dtw) {
                        base_best_case += fr_rt + hd_rt/1000 + dtw_rt/1000;
                        base_worst_case += fr_rt + hd_rt/1000 + dtw_rt/1000;
                    }
                    else{
                        (fr) ? base_worst_case += fr_rt : curr_best_case = fr_rt, curr_worst_case = fr_rt;
                        (hd) ? base_worst_case += hd_rt/1000 : curr_best_case = std::min(curr_best_case, hd_rt/1000), curr_worst_case = std::max(curr_worst_case, hd_rt/1000);
                        (dtw) ? base_worst_case += dtw_rt/1000 : curr_best_case = std::min(curr_best_case, dtw_rt/1000), curr_worst_case = std::max(curr_worst_case, dtw_rt/1000);
                        base_best_case += curr_best_case;
                        base_worst_case += curr_worst_case;
                    }

                    break;
                case OR:

                    if((!fr) && (!hd) && (!dtw)){
                        base_best_case += fr_rt + hd_rt/1000 + dtw_rt/1000;
                        base_worst_case += fr_rt + hd_rt/1000 + dtw_rt/1000;
                    }
                    else{
                        (fr) ? curr_best_case = fr_rt, curr_worst_case = fr_rt : base_worst_case += fr_rt;
                        (hd) ? curr_best_case = std::min(curr_best_case, hd_rt/1000), curr_worst_case = std::max(curr_worst_case, hd_rt/1000) : base_worst_case += hd_rt/1000;
                        (dtw) ? curr_best_case = std::min(curr_best_case, dtw_rt/1000), curr_worst_case = std::max(curr_worst_case, dtw_rt/1000) : base_worst_case += dtw_rt/1000;
                        base_best_case += curr_best_case;
                        base_worst_case += curr_worst_case;
                    }

                    break;
                case MAJORITY:

                    if(int(fr) + int(hd) + int(dtw) == 3 || int(fr) + int(hd) + int(dtw) == 0){

                        double rts[3] = {fr_rt, hd_rt/1000, dtw_rt/1000};
                        std::sort(rts, rts + 3);

                        base_worst_case += rts[1] + rts[2];                         // worst = two maximal runtimes
                        base_best_case += rts[0] + rts[1];                         // best = two minimal runtimes

                    }else if(int(fr) + int(hd) + int(dtw) == 2) {

                        base_worst_case += fr_rt + hd_rt/1000 + dtw_rt/1000;                // worst = rt of all three

                        (fr) ?      ((hd) ? (base_best_case += fr_rt + hd_rt/1000) : (base_best_case += fr_rt + dtw_rt/1000) )
                                  : base_best_case += hd_rt/1000 + dtw_rt/1000;           // best = rt of the two positive results


                    }else if(int(fr) + int(hd) + int(dtw) == 1){

                        base_worst_case += fr_rt + hd_rt/1000 + dtw_rt/1000;                // worst = rt of all three

                        (!fr) ?      ((!hd) ? (base_best_case += fr_rt + hd_rt/1000) : (base_best_case += fr_rt + dtw_rt/1000))
                             : base_best_case += hd_rt/1000 + dtw_rt/1000;                 // best = rt of the two negative results


                    }
                    break;
            }

            compared_trajectories++;

        }
    }

    // Iterate over trajectories - Combined Computation
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {

            std::tuple<double, double, double> runtimes;

            runtimes = timed_fr_hd_dtw_similarity_join(**i_it, **j_it, epsilons[0], epsilons[2], epsilons[1],combiner,false);
            comb_fr += get<0>(runtimes);
            comb_hd += get<1>(runtimes);
            comb_dtw += get<2>(runtimes);

            double rt = get<0>(runtimes) + get<1>(runtimes)/1000 + get<2>(runtimes)/1000;
            comb_rt += rt;
        }
    }

    std::cout << "Done!" << std::endl;
    std::cout << "Trajectory pairs checked: " << compared_trajectories << std::endl;
    std::cout << std::endl;
    std::cout << "Results:" << std::endl << std::endl;
    std::cout << COUT_TAB << "Best Indep. RT." << COUT_TAB << "Worst Indep. RT" << COUT_TAB << "Best Ind. RT. Avg." << COUT_TAB << "Worst Ind. RT. Avg." << COUT_TAB << "Comb. RT." << COUT_TAB << "Comb. RT. Avg." << std::endl;
    std::cout << COUT_TAB << base_best_case << COUT_TAB << base_worst_case << COUT_TAB << base_best_case/compared_trajectories << COUT_TAB << base_worst_case/compared_trajectories << COUT_TAB << comb_rt << COUT_TAB << comb_rt/compared_trajectories << std::endl;
    std::cout << std::endl << std::endl;
    std::cout << COUT_TAB << "Indep. Fr" << COUT_TAB << "Indep. Hd." << COUT_TAB << "Indep. DTW" << COUT_TAB << "Comb. Fr" << COUT_TAB << "Comb. Hd." << COUT_TAB << "Comb. DTW" << std::endl;
    std::cout << COUT_TAB << base_fr << COUT_TAB << base_hd/1000 << COUT_TAB << base_dtw/1000 << COUT_TAB << comb_fr << COUT_TAB << comb_hd/1000 << COUT_TAB << comb_dtw/1000 << std::endl;
    std::cout << std::endl << std::endl;
}

template<class T, class DataLoader>
void CombinedMeasuresExperiments<T, DataLoader>::fs_exploral(unsigned int trajectory_cardinality, std::array<float, 3> epsilons) {

    std::cout << "Running Experiment: Effect of Explored Free Space" << std::endl;

    // Load Data
    std::cout << "Loading Dataset...";
    std::shared_ptr<std::shared_ptr<Trajectory<T>>[]> trajectories;
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Setup:" << std::endl;
    std::cout << "n = " << trajectory_cardinality << "; eps_fr = " << epsilons[0] << "; eps_dtw = " << epsilons[1] << "; eps_hd = " << epsilons[2] << std::endl;

    std::cout << "Running Tests...";

    // Setup Variables
    double fr_exploration = 0, comb_rt = 0, fr_rt = 0, dtw_rt = 0, hd_rt = 0;
    int nr_trajectories = 0;

    // Iterate over trajectories
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {

            auto runtimes = timed_exploral_fr_hd_dtw_algorithm(**i_it, **j_it, epsilons[0], epsilons[2], epsilons[1], true);

            if((**i_it).length() != 0 && (**j_it).length() != 0){
                fr_exploration += get<0>(runtimes) / float((**i_it).length() * (**j_it).length());
            }

            fr_rt += get<1>(runtimes);
            hd_rt += get<2>(runtimes);
            dtw_rt += get<3>(runtimes);
            comb_rt += get<1>(runtimes) + get<2>(runtimes)/1000 + get<3>(runtimes)/1000;

            nr_trajectories++;
        }
    }

    std::cout << "Done!" << std::endl;
    std::cout << "Tested Trajectories: " << nr_trajectories << std::endl;
    std::cout << std::endl;
    std::cout << "Results:" << std::endl << std::endl;
    std::cout << COUT_TAB << "Explored FS. Avg." << COUT_TAB << "Total RT." << COUT_TAB << "Fr. RT." << COUT_TAB << "Hd. RT." << COUT_TAB << "DTW. RT." << std::endl;
    std::cout << COUT_TAB << fr_exploration/float(nr_trajectories) << COUT_TAB << comb_rt << COUT_TAB << fr_rt << COUT_TAB << hd_rt/1000 << COUT_TAB << dtw_rt/1000 << std::endl;
    std::cout << COUT_TAB << "Explored FS. Avg." << COUT_TAB << "Total RT. Avg." << COUT_TAB << "Fr. RT. Avg." << COUT_TAB << "Hd. RT. Avg." << COUT_TAB << "DTW. RT. Avg." << std::endl;
    std::cout << COUT_TAB << fr_exploration/float(nr_trajectories) << COUT_TAB << comb_rt/float(nr_trajectories) << COUT_TAB << fr_rt/float(nr_trajectories) << COUT_TAB << hd_rt/float(nr_trajectories * 1000) << COUT_TAB << dtw_rt/float(nr_trajectories * 1000) << std::endl;


}

template<class T, class DataLoader>
void CombinedMeasuresExperiments<T, DataLoader>::selectivity(unsigned int trajectory_cardinality,std::array<float, 3> epsilons) {

    std::cout << "Running Experiment: \"Selectivity\" of Combined Similarity Measure" << std::endl;

    // Load Data
    std::cout << "Loading Dataset...";
    std::shared_ptr<std::shared_ptr<Trajectory<T>>[]> trajectories;
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Setup:" << std::endl;
    std::cout << "n = " << trajectory_cardinality << "; eps_fr = " << epsilons[0] << "; eps_dtw = " << epsilons[1] << "; eps_hd = " << epsilons[2] << std::endl;

    std::cout << "Running Tests...";

    // Setup Variables
    int fr_pos = 0, hd_pos = 0, dtw_pos = 0, and_pos = 0, or_pos = 0, majority_pos = 0, nr_trajectories = 0;

    // Iterate over trajectories
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {

                bool fr,hd,dtw;

                auto cache = std::make_shared<DistCache<T>>(*i_it, *j_it);
                fr = frechet_LT_Epsilon(**i_it, **j_it, epsilons[0], cache);
                dtw = dtw_LT_epsilon(**i_it, **j_it, epsilons[1], cache);
                hd = hd_LT_epsilon(**i_it, **j_it, epsilons[2], cache);

                if(fr) {fr_pos++;}
                if(dtw) {dtw_pos++;}
                if(hd) {hd_pos++;}

                if(fr && dtw && hd) {and_pos++;}
                if(fr || dtw || hd) {or_pos++;}
                if(int(fr) + int(dtw) + int(hd) >= 2) {majority_pos++;}

                nr_trajectories++;

        }
    }

    std::cout << "Done!" << std::endl;
    std::cout << "Tested Trajectories: " << nr_trajectories << std::endl;
    std::cout << std::endl;
    std::cout << "Results:" << std::endl << std::endl;
    std::cout << COUT_TAB << "Fr. Sel." << COUT_TAB << "Hd. Sel." << COUT_TAB << "DTW. Sel." << COUT_TAB << "And Sel." << COUT_TAB << "Or Sel." << COUT_TAB << "Maj. Sel." << std::endl;
    std::cout << COUT_TAB << float(fr_pos)/nr_trajectories << COUT_TAB << float(hd_pos)/nr_trajectories << COUT_TAB << float(dtw_pos)/nr_trajectories << COUT_TAB << float(and_pos)/nr_trajectories << COUT_TAB << float(or_pos)/nr_trajectories << COUT_TAB << float(majority_pos)/nr_trajectories << std::endl << std::endl;

}

template<class T, class DataLoader>
void CombinedMeasuresExperiments<T, DataLoader>::pure_runtime_with_dfr(unsigned int trajectory_cardinality, std::array<float, 3> epsilons) {

    std::cout << "Running Experiment: Combined vs. Independent Computation - Pure Runtime" << std::endl;

    // Load Data
    std::cout << "Loading Dataset...";
    std::shared_ptr<std::shared_ptr<Trajectory<T>>[]> trajectories;
    unsigned int n = -1;
    std::throw_with_nested("Dead Code. To use this function, please implement a call to an according dataset or adapt the code.");
    std::cout << "Done!" << std::endl;

    std::cout << "Setup:" << std::endl;
    std::cout << "n = " << trajectory_cardinality << "; eps_fr = " << epsilons[0] << "; eps_dtw = " << epsilons[1] << "; eps_hd = " << epsilons[2] << std::endl;

    std::cout << "Running Tests...";

    // Setup Runtime Variables
    double base_runtime = 0, comb_runtime = 0, base_rt_best = INFINITY, comb_rt_best = INFINITY, base_rt_worst = 0, comb_rt_worst = 0;
    double base_fr = 0, base_hd = 0, base_dtw = 0, comb_fr = 0, comb_hd = 0, comb_dtw = 0;
    int compared_trajectories = 0;

    // Iterate over trajectories - Base Computation
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {

            double fr_rt = 0, hd_rt = 0, dtw_rt = 0;
            auto cache = std::make_shared<DistCache<T>>(*i_it, *j_it);

            auto start = std::chrono::high_resolution_clock::now();
            dFr_LT_epsilon(**i_it, **j_it, epsilons[0], cache);
            auto end = std::chrono::high_resolution_clock::now();
            fr_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            base_fr += fr_rt;

            start = std::chrono::high_resolution_clock::now();
            hd_LT_epsilon(**i_it,**j_it,epsilons[2], cache);
            end = std::chrono::high_resolution_clock::now();
            hd_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            base_hd += hd_rt;

            start = std::chrono::high_resolution_clock::now();
            dtw_LT_epsilon(**i_it,**j_it,epsilons[1], cache);
            end = std::chrono::high_resolution_clock::now();
            dtw_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            base_dtw += dtw_rt;


            double rt = (fr_rt + hd_rt + dtw_rt);
            base_runtime += rt;

            if(rt < base_rt_best) {base_rt_best = rt;}
            if(rt > base_rt_worst) {base_rt_worst = rt;}

            compared_trajectories++;

        }
    }

    // Iterate over trajectories - Combined Computation
    for (auto i_it = trajectories.get(); i_it < trajectories.get() + n; i_it++) {
        for (auto j_it = i_it + 1; j_it < trajectories.get() + n; j_it++) {


            auto start = std::chrono::high_resolution_clock::now();
            dfr_hd_dtw_algorithm(**i_it, **j_it, epsilons[0], epsilons[2], epsilons[1], true);
            auto end = std::chrono::high_resolution_clock::now();


            double rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            comb_runtime += rt;

            if(rt < comb_rt_best) {comb_rt_best = rt;}
            if(rt > comb_rt_worst) {comb_rt_worst = rt;}
        }
    }

    std::cout << "Done!" << std::endl;
    std::cout << "Trajectory pairs checked: " << compared_trajectories << std::endl;
    std::cout << std::endl;
    std::cout << "Results:" << std::endl << std::endl;

    std::cout << COUT_TAB << "RT. Indep." << COUT_TAB << "RT. Indep. Avg." << COUT_TAB << "RT. Indep. Best" << COUT_TAB << "RT. Indep. Worst" << COUT_TAB << "RT. Comb." << COUT_TAB << "RT. Comb. Avg." << COUT_TAB << "RT. Comb. Best" << COUT_TAB << "RT. Comb. Worst." << std::endl;
    std::cout << COUT_TAB << base_runtime/1000 << COUT_TAB << base_runtime/(1000 * compared_trajectories) << COUT_TAB << base_rt_best/1000 << COUT_TAB << base_rt_worst/1000 << COUT_TAB << comb_runtime/1000 << COUT_TAB << comb_runtime/(1000 * compared_trajectories) << COUT_TAB << comb_rt_best/1000 << COUT_TAB << comb_rt_worst/1000 << std::endl;
    std::cout << std::endl;

}