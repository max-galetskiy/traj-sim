#pragma once
#include "../TrajModels/RoadNetworkPoint.hpp"
#include "../TrajModels/Trajectory.hpp"

class DataLoaderRoadNetwork {

    public:

    inline static std::shared_ptr<std::shared_ptr<Trajectory<RoadNetworkPoint>>[]> trajectories;
    inline static unsigned int trajectories_size = 0;

    static void readFile(std::string graph_filepath, std::string trajectory_filepath, std::string phl_infilepath, std::string phl_outfilepath, unsigned int traj_cardinality, bool graph_header_present, bool traj_header_present, unsigned int source_column, unsigned int target_column, unsigned int distance_column, unsigned int traj_column);

    static Graph readGraph(std::string filepath, bool header_present, unsigned int source_column, unsigned int target_column, unsigned int distance_column);
    static std::shared_ptr<PrunedHighwayLabeling> readPHL(std::string input_filepath, std::string output_filepath);

    static void readTrajectories(std::string filepath, unsigned int traj_cardinality, bool header_present, unsigned int traj_column);

};

