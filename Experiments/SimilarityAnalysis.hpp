#pragma once
#include "../DataLoading/DataLoaderEuclidean.hpp"
#include "../DataLoading/DataLoaderRoadNetwork.hpp"

typedef std::shared_ptr<std::shared_ptr<Trajectory<EuclideanPoint2D>>[]> euclidean_trajectories;
typedef std::shared_ptr<std::shared_ptr<Trajectory<RoadNetworkPoint>>[]> roadnet_trajectories;

class SimilarityAnalysis {

    private:
        unsigned int traj_cardinality = 0;
        euclidean_trajectories eucl_trajs;
        roadnet_trajectories roadnet_trajs;

    public:

        explicit SimilarityAnalysis(unsigned int traj_cardinality);
        void compare_to_jaccard(float tau);
        static float jaccard_similarity(Trajectory<RoadNetworkPoint>& t1, Trajectory<RoadNetworkPoint>& t2);

};



