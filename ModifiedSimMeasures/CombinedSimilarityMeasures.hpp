#pragma once

#include <chrono>

#include "ModifiedFrechet.hpp"
#include "ModifiedHausdorff.hpp"
#include "ModifiedDTW.hpp"

enum COMB_MEASURE {AND, OR, MAJORITY};


template<class T>
std::tuple<bool,bool,bool> frechet_hausdorff_dtw_algorithm(Trajectory<T>& t1, Trajectory<T>& t2, float fr_eps, float hd_eps, float dtw_eps, bool use_cache = false){

    bool fr_LT_eps, hd_LT_eps, dtw_LT_eps;

    auto hd_ind = HausdorffIndexManager(t1.length(), t2.length());
    auto dtw_paths = DTWPathManager(t1.length(), t2.length());

    auto cache = (use_cache) ? std::make_shared<DistCache<T>>(std::make_shared<Trajectory<T>>(t1),std::make_shared<Trajectory<T>>(t2)) : nullptr;

    fr_LT_eps = modified_frechet_LT_epsilon(t1,t2,fr_eps, hd_eps, hd_ind, dtw_paths, cache);

    if(hd_ind.decided()){
        hd_LT_eps = hd_ind.lt_eps();
    }else{
        hd_LT_eps = modified_hd_LT_epsilon(t1,t2,hd_eps, hd_ind, cache);
    }

    if(auto dtw_bound = dtw_paths.get_bound(t1.length() - 1,t2.length() - 1); dtw_bound.second != NONE){
        if(dtw_bound.second == GREATER && dtw_bound.first > dtw_eps){
            dtw_LT_eps = false;
        }
        else if(dtw_bound.second == LESS_EQ && dtw_bound.first <= dtw_eps){
            dtw_LT_eps = true;
        }
        else{
            dtw_LT_eps = modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
        }

    }else{
        dtw_LT_eps = modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
    }

    return {fr_LT_eps,hd_LT_eps,dtw_LT_eps};
}


template<class T>
bool fr_hd_dtw_similarity_join(Trajectory<T>& t1, Trajectory<T>& t2, float fr_eps, float hd_eps, float dtw_eps, COMB_MEASURE combiner = AND, bool use_cache = false){
    bool fr_LT_eps, hd_LT_eps, dtw_LT_eps;

    auto hd_ind = HausdorffIndexManager(t1.length(), t2.length());
    auto dtw_paths = DTWPathManager(t1.length(), t2.length());

    auto cache = (use_cache) ? std::make_shared<DistCache<T>>(std::make_shared<Trajectory<T>>(t1),std::make_shared<Trajectory<T>>(t2)) : nullptr;

    fr_LT_eps = modified_frechet_LT_epsilon(t1,t2,fr_eps, hd_eps, hd_ind, dtw_paths, cache);

    // Early Breaks
    if(combiner == AND && !fr_LT_eps){ return false;}
    else if(combiner == OR && fr_LT_eps){return true;}

    if(hd_ind.decided()){
        hd_LT_eps = hd_ind.lt_eps();
    }else{
        hd_LT_eps = modified_hd_LT_epsilon(t1,t2,hd_eps, hd_ind, cache);
    }

    // Early Breaks
    if(combiner == AND && !hd_LT_eps){return false;}
    else if(combiner == OR && hd_LT_eps){return true;}
    else if(combiner == MAJORITY && (fr_LT_eps == hd_LT_eps)){return fr_LT_eps;}

    if(auto dtw_bound = dtw_paths.get_bound(t1.length() - 1,t2.length() - 1); dtw_bound.second != NONE){
        if(dtw_bound.second == GREATER && dtw_bound.first > dtw_eps){
            dtw_LT_eps = false;
        }
        else if(dtw_bound.second == LESS_EQ && dtw_bound.first <= dtw_eps){
            dtw_LT_eps = true;
        }
        else{
            dtw_LT_eps = modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
        }

    }else{
        dtw_LT_eps = modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
    }

    if(combiner == AND || combiner == OR){return dtw_LT_eps;}
    else{
        return (int(fr_LT_eps) + int(dtw_LT_eps) + int(hd_LT_eps)) >= 2;
    }

}


template<class T>
std::tuple<bool,bool,bool> dfr_hd_dtw_algorithm(Trajectory<T>& t1, Trajectory<T>& t2, float fr_eps, float hd_eps, float dtw_eps, bool use_cache = false){
    bool fr_LT_eps, hd_LT_eps, dtw_LT_eps;

    auto hd_ind = HausdorffIndexManager(t1.length(), t2.length());
    auto dtw_paths = DTWPathManager(t1.length(), t2.length());

    auto cache = (use_cache) ? std::make_shared<DistCache<T>>(std::make_shared<Trajectory<T>>(t1),std::make_shared<Trajectory<T>>(t2)) : nullptr;

    fr_LT_eps = modified_discrete_frechet_LT_epsilon(t1,t2,fr_eps, hd_eps, hd_ind, dtw_paths, cache);

    if(hd_ind.decided()){
        hd_LT_eps = hd_ind.lt_eps();
    }else{
        hd_LT_eps = modified_hd_LT_epsilon(t1,t2,hd_eps, hd_ind, cache);
    }

    if(auto dtw_bound = dtw_paths.get_bound(t1.length() - 1,t2.length() - 1); dtw_bound.second != NONE){
        if(dtw_bound.second == GREATER && dtw_bound.first > dtw_eps){
            dtw_LT_eps = false;
        }
        else if(dtw_bound.second == LESS_EQ && dtw_bound.first <= dtw_eps){
            dtw_LT_eps = true;
        }
        else{
            dtw_LT_eps = modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
        }

    }else{
        dtw_LT_eps = modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
    }

    return {fr_LT_eps,hd_LT_eps,dtw_LT_eps};
};