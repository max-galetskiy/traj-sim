#pragma once
#include <vector>
#include <stdexcept>
#include "../TrajModels/Trajectory.hpp"

/*
    Class for a cache of point-wise distances (potentially aiding the evaluation of a combined similarity measure)


	Template class Point needs to support the following functions: 
		- static float Point::dist(Point& a, Point& b)
*/ 
template<class Point>
class DistCache {
    private:

        std::unique_ptr<float[]> dists;
        std::unique_ptr<bool[]> dists_cached;

        std::shared_ptr<Trajectory<Point>> t1;
        std::shared_ptr<Trajectory<Point>> t2;
        unsigned int n;
        unsigned int m;

    public:

        // Optionally fills the entire distance cache immediately
        DistCache(std::shared_ptr<Trajectory<Point>> t1, std::shared_ptr<Trajectory<Point>> t2, bool fill = false) : t1(t1), t2(t2) {
            n  = (*t1).length();
            m  = (*t2).length();

            // Immediately fill all
            if(fill){
                dists = std::move(cross_product(*t1, *t2));
                dists_cached = std::unique_ptr<bool[]>(new bool[n * m]);
                for(int i = 0; i < n*m; i++){
                    dists_cached[i] = true;
                }

            // Initialize empty
            }
            else{
                dists_cached = std::make_unique<bool[]>(n * m);         // Fills the array with 0s! 
                dists = std::unique_ptr<float[]>(new float[n * m]);          // Does not fill the array!
            }

        }

        // Computes/Retrieves distance between points i and j
        float cacheDistance(unsigned int i, unsigned int j){

            if(i >= n || j >= m){
                std::throw_with_nested(std::invalid_argument("Distance Cache Index out of Bounds!"));
            }

            if(dists_cached[i * m + j]){
                return dists[i * m + j];
            }
            else{
                float d = Point::dist((*t1).at(i), (*t2).at(j));
                dists[i * m + j] = d;
                dists_cached[i * m + j] = true;
                return d;
            }

        }

        // Computes all point-wise distances of two trajectories
        static std::unique_ptr<float[]> cross_product(Trajectory<Point>& t1, Trajectory<Point>& t2){

            const unsigned int n = t1.length();
            const unsigned int m = t2.length();

            std::unique_ptr<float[]> d(new float[n * m]);

            for (unsigned int i = 0; i < n; i++){

                for(unsigned int j = 0; j < m; j++){

                    d[i * m + j] = Point::dist(t1.at(i), t2.at(j));

                }
            }

            return d;

        }
};