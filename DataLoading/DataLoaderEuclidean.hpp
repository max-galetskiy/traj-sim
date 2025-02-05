#pragma once
#include <string>
#include <memory>

#include "../TrajModels/Trajectory.hpp"
#include "../TrajModels/EuclideanPoint2D.hpp"

class DataLoaderEuclidean {

public:

    inline static std::shared_ptr<std::shared_ptr<Trajectory<EuclideanPoint2D>>[]> trajectories;
    inline static unsigned int trajectories_size = 0;

    static void readFile(const std::string& filepath, unsigned int traj_cardinality, bool header_present, unsigned int traj_column);

};

