#include "SimilarityAnalysis.hpp"
#include <fstream>
#include "../SimMeasures/DiscreteFrechet.hpp"
#include "../SimMeasures/Hausdorff.hpp"
#include "../SimMeasures/DTW.hpp"
#include "../SimMeasures/EDR.hpp"
#include "../SimMeasures/LCSS.hpp"
#include "../Utils/Utilities.hpp"

SimilarityAnalysis::SimilarityAnalysis(unsigned int traj_cardinality) : traj_cardinality(traj_cardinality) {

    DataLoaderEuclidean::readFile("../Datasets/euclidean_sample.csv",4,true,1);
    DataLoaderRoadNetwork::readFile("../Datasets/rn_graph_sample.csv","../Datasets/rn_sample.csv","../Datasets/rn_graph_sample_phl.txt","../Datasets/rn_graph_sample_labels.txt",3,true,true,0,1,2,1);
    eucl_trajs = DataLoaderEuclidean::trajectories;
    roadnet_trajs = DataLoaderRoadNetwork::trajectories;
}

/*
    For each trajectory pair we compare the jaccard similarity of the roadnet graph to our similarity measures
 */
void SimilarityAnalysis::compare_to_jaccard(float tau) {

    std::ofstream log("../Experiments/sim_analysis_data/sim_values_with_jaccard.csv",std::fstream::out);

    if(!log){
        std::throw_with_nested(std::runtime_error("Could not create sim_values_with_jaccard.csv!"));
        return;
    }

    log << "i;j;jaccard;dFr;Hd;DTW;EDR;LCSS" << std::endl;

    for(int i = 0; i < traj_cardinality; i++){
       for(int j = i + 1; j < traj_cardinality; j++){

           auto t1 = eucl_trajs.get()[i];
           auto t2 = eucl_trajs.get()[j];

           auto t1_matched = roadnet_trajs.get()[i];
           auto t2_matched = roadnet_trajs.get()[j];

           log << i << ";" << j << ";";

           float jaccard = SimilarityAnalysis::jaccard_similarity(*t1_matched,*t2_matched);
           log << jaccard << ";";

           log  << discrete_frechet(*t1,*t2) << ";"
                << hausdorff(*t1,*t2) << ";"
                << dynamic_time_warping(*t1,*t2) << ";"
                << edit_distance_real_sequence(*t1,*t2,tau) << ";"
                << longest_common_subsequence(*t1,*t2,tau)  << std::endl;

       }
    }

    log.close();
}

/*
 *      Computes the jaccard similarity between two RoadNet trajectories defined as:
 *       |E1 n E2| / |E1 u E2|   or equivalently: |E1 n E2| / (|E1| + |E2| - |E1 n E2|)
 *      This code utilizes the property that sequential nodes in a road network trajectory are always edges in the road network graph
 */
float SimilarityAnalysis::jaccard_similarity(Trajectory<RoadNetworkPoint>& t1, Trajectory<RoadNetworkPoint>& t2) {

    std::unordered_set<std::pair<int,int>,pair_hash> seen_edges_t1;
    std::unordered_set<std::pair<int,int>,pair_hash> seen_edges_t2;
    unsigned int t1_edges = 0;
    unsigned int t2_edges = 0;
    unsigned int intersect_edges = 0;

    if(t1.length() < 2 || t2.length() < 2){     // No edges exist
        //std::cout << t1.length() << " : " << t2.length() << std::endl;
        return -1;
    }

    // NOTE: WE IGNORE THE FIRST AND LAST EDGES
    for(int i = 1; i < t1.length() - 2; i++){

        auto edge = std::pair<int,int>(t1.at(i).node_id,t1.at(i+1).node_id);

        if(seen_edges_t1.find(edge) == seen_edges_t1.end()){                                      // We want to avoid counting the same edge twice
            seen_edges_t1.insert(edge);
            t1_edges++;
        }
    }

    // NOTE: We ignore first and last edges
    for(int i = 1; i < t2.length() - 2; i++){

        auto edge = std::pair<int,int>(t2.at(i).node_id,t2.at(i+1).node_id);

        if(seen_edges_t2.find(edge) == seen_edges_t2.end()){
            seen_edges_t2.insert(edge);
            t2_edges++;

            if(seen_edges_t1.find(edge) != seen_edges_t1.end()){
                intersect_edges++;
            }

        }

    }

    return (intersect_edges * 1.f / float(t1_edges + t2_edges - intersect_edges));

}
