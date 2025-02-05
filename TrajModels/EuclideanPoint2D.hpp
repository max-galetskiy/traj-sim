#pragma once
#include <string>
#include <cmath>

#include "FreeSpaceCell.hpp"
#include "../DistCache/DistCache.hpp"

/*
    Class representing GPS points in a 2D Euclidean space
    Contains multiple 2D-geometry-related functions as well
 */
class EuclideanPoint2D {
private:

    static fpair circleLineIntersection(EuclideanPoint2D& center, EuclideanPoint2D& l1, EuclideanPoint2D& l2, float r, int i);

public:

    static const int fs_intervals_per_side = 1;

	float x;
	float y;

	EuclideanPoint2D(float x, float y) : x(x), y(y) {};
    EuclideanPoint2D(std::pair<float,float> longLat, bool meters);

	std::string to_string();

	static float dist(EuclideanPoint2D& a, EuclideanPoint2D& b);
    static float dist_to_edge(EuclideanPoint2D& p, EuclideanPoint2D& e1, EuclideanPoint2D &e2);
    static float dist_of_bisector(EuclideanPoint2D& p1, EuclideanPoint2D& p2, EuclideanPoint2D& e1, EuclideanPoint2D& e2);

    static FreeSpaceCell createFreeSpace(EuclideanPoint2D& p1, EuclideanPoint2D& q1, EuclideanPoint2D& p2, EuclideanPoint2D& q2, float epsilon, int i, int j, std::shared_ptr<DistCache<EuclideanPoint2D>> cache = nullptr);
    static void propagateFreeSpace(EuclideanPoint2D& p1, EuclideanPoint2D& q1, EuclideanPoint2D& p2, EuclideanPoint2D& q2, FreeSpaceCell& cell, float epsilon, std::shared_ptr<DistCache<EuclideanPoint2D>> cache = nullptr);
    static void computeFreeSpaceSide(EuclideanPoint2D &p1, EuclideanPoint2D &q1, EuclideanPoint2D &p2, EuclideanPoint2D &q2, FreeSpaceCell &cell, float epsilon, side s, std::shared_ptr<DistCache<EuclideanPoint2D>> cache = nullptr);
};
