#pragma once

#include <variant>

#include "../TrajModels/Trajectory.hpp"
#include "../DistCache/DistCache.hpp"

/*
    This class mainly deals with the different parameter types of the different similarity measures
    (Without it, code for baselines becomes quite messy and redundant)
*/
template<class T>
class AbstractMeasure {
private:

    typedef float (*normal_measure)(Trajectory<T> &,Trajectory<T> &,std::shared_ptr<DistCache<T>>);
    typedef float (*ed_measure)(Trajectory<T> &, Trajectory<T> &, float, std::shared_ptr<DistCache<T>>);
    typedef float (*erp_measure)(Trajectory<T> &, Trajectory<T> &, T, float, std::shared_ptr<DistCache<T>>);
    typedef float (*rlcss_measure)(Trajectory<T> &, Trajectory<T> &, float, int, std::shared_ptr<DistCache<T>>);

    typedef bool (*normal_decider)(Trajectory<T> &,Trajectory<T> &, float, std::shared_ptr<DistCache<T>>);
    typedef bool (*ed_decider)(Trajectory<T> &, Trajectory<T> &, float, float, std::shared_ptr<DistCache<T>>);
    typedef bool (*erp_decider)(Trajectory<T> &, Trajectory<T> &, T, float, float, std::shared_ptr<DistCache<T>>);
    typedef bool (*rlcss_decider)(Trajectory<T> &, Trajectory<T> &, float, int, float, std::shared_ptr<DistCache<T>>);

    static inline float match_threshold = -1;
    static inline int index_threshold = 0;
    static inline std::shared_ptr<T> ref_point_ptr = nullptr;

    static inline std::variant<normal_measure,ed_measure,erp_measure,rlcss_measure> dist_function = {};
    static inline std::variant<normal_decider, ed_decider , erp_decider , rlcss_decider> dist_decider = {};

public:

    /*
            Computation
     */
    static void setMeasure(normal_measure f){
        dist_function = f;
        calc_dist = std::get<normal_measure>(dist_function);
    }

    static void setMeasure(ed_measure f, float tau){
        dist_function = f;
        match_threshold = tau;
        calc_dist = [](Trajectory<T>& t1, Trajectory<T>& t2, std::shared_ptr<DistCache<T>> c){ return (std::get<ed_measure>(dist_function))(t1,t2,match_threshold, c); };
    }

    static void setMeasure(erp_measure f, float tau, T g){
        dist_function = f;
        match_threshold = tau;
        ref_point_ptr = std::make_shared<T>(g);
        calc_dist = [](Trajectory<T>& t1, Trajectory<T>& t2, std::shared_ptr<DistCache<T>> c){ return (std::get<erp_measure>(dist_function))(t1,t2,*ref_point_ptr,match_threshold, c); };
    }

    static void setMeasure(rlcss_measure f, float tau, int delta){
        dist_function = f;
        match_threshold = tau;
        index_threshold = delta;
        calc_dist = [](Trajectory<T>& t1, Trajectory<T>& t2, std::shared_ptr<DistCache<T>> c){ return (std::get<rlcss_measure >(dist_function))(t1,t2, match_threshold, index_threshold, c); };
    }

    static inline float(*calc_dist)(Trajectory<T>&, Trajectory<T>&, std::shared_ptr<DistCache<T>>) = 0;

    /*
            Decision Problem
    */
    static void setDecider(normal_decider f){
        dist_decider = f;
        decide_dist = std::get<normal_decider>(dist_decider);
    }

    static void setDecider(ed_decider f, float tau){
        dist_decider = f;
        match_threshold = tau;
        decide_dist = [](Trajectory<T>& t1, Trajectory<T>& t2, float e, std::shared_ptr<DistCache<T>> c){ return (std::get<ed_decider>(dist_decider))(t1,t2,match_threshold, e, c); };
    }

    static void setDecider(erp_decider f, float tau, T g){
        dist_decider = f;
        match_threshold = tau;
        ref_point_ptr = std::make_shared<T>(g);
        decide_dist = [](Trajectory<T>& t1, Trajectory<T>& t2, float e, std::shared_ptr<DistCache<T>> c){ return (std::get<erp_decider>(dist_decider))(t1,t2,*ref_point_ptr,match_threshold, e, c); };
    }

    static void setDecider(rlcss_decider f, float tau, int delta){
        dist_decider = f;
        match_threshold = tau;
        index_threshold = delta;
        decide_dist = [](Trajectory<T>& t1, Trajectory<T>& t2, float e, std::shared_ptr<DistCache<T>> c){ return (std::get<rlcss_decider >(dist_decider))(t1,t2, match_threshold, index_threshold, e, c); };
    }

    static inline bool(*decide_dist)(Trajectory<T>&, Trajectory<T>&, float, std::shared_ptr<DistCache<T>>) = 0;
};

