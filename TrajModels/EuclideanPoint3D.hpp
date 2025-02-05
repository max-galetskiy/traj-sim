#pragma once
#include <string>

#include "../DistCache/DistCache.hpp"
#include "FreeSpaceCell.hpp"

/*
    Adapts the two-dimensional model to a 3D Euclidean space
 */
class EuclideanPoint3D {
private:

    static fpair circleLineIntersection(EuclideanPoint3D& center, EuclideanPoint3D& l1, EuclideanPoint3D& l2, float r, int i);

public:

    static const int fs_intervals_per_side = 1;

    float x;
    float y;
    float z;

    EuclideanPoint3D(float x, float y, float z) : x(x), y(y), z(z) {};
    EuclideanPoint3D(std::pair<float,float> longLat, bool meters = false);

    std::string to_string();

    static float dist(EuclideanPoint3D& a, EuclideanPoint3D& b);
    static float dist_to_edge(EuclideanPoint3D& p, EuclideanPoint3D& e1, EuclideanPoint3D &e2);
    static float dist_of_bisector(EuclideanPoint3D& p1, EuclideanPoint3D& p2, EuclideanPoint3D& e1, EuclideanPoint3D& e2);

    static FreeSpaceCell createFreeSpace(EuclideanPoint3D& p1, EuclideanPoint3D& q1, EuclideanPoint3D& p2, EuclideanPoint3D& q2, float epsilon, int i, int j, std::shared_ptr<DistCache<EuclideanPoint3D>> cache = NULL);
    static void propagateFreeSpace(EuclideanPoint3D& p1, EuclideanPoint3D& q1, EuclideanPoint3D& p2, EuclideanPoint3D& q2, FreeSpaceCell& cell, float epsilon, std::shared_ptr<DistCache<EuclideanPoint3D>> cache = NULL);
    static void computeFreeSpaceSide(EuclideanPoint3D &p1, EuclideanPoint3D &q1, EuclideanPoint3D &p2, EuclideanPoint3D &q2, FreeSpaceCell &cell, float epsilon, side s, std::shared_ptr<DistCache<EuclideanPoint3D>> cache = NULL);

};


