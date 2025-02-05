#pragma once 
#include <deque>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <iostream>

#include "../TrajModels/Trajectory.hpp"
#include "../TrajModels/FreeSpaceCell.hpp"
#include "../DistCache/DistCache.hpp"

#include "../Utils/Bitmask.hpp"

/*
    Note that Frechet here refers to what some call the continuous Frechet distance (as opposed to the discrete one!)

    class T needs to support: 
        - static float T::dist(T& p1, T& p2)
        - static float T::dist_to_edge(T& p, T& e1, T& e2)
        - static float T::dist_of_bisector(T& p1, T& p2, T& e1, T& e2)
        - static FreeSpaceCell T::createFreeSpace(T& p1, T& q1, T& p2, T& q2, float epsilon, int i, int j, std::shared_ptr<DistCache<T>> cache = NULL)
                    // Creates the FreeSpace Diagram normally
        - static void T::propagateFreeSpace(T& p1, T& q1, T& p2, T& q2, FreeSpace& cell, float epsilon, std::shared_ptr<DistCache<T>> cache = NULL)
                    // Adjusts the right and top intervals to be valid for the given Free Space cell
        - static void T::computeFreeSpaceSide(T& p1, T& q1, T& p2, T& q2, FreeSpace& cell, float epsilon, Side s, std::shared_ptr<DistCache<T>> cache = NULL)
                    // Computes the given side for the given Free Space cell
        - static const int fs_intervals_per_side
                    // Number of intervals on a side of the Free Space Diagram
*/
template<class T>
bool frechet_LT_Epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float epsilon, std::shared_ptr<DistCache<T>> cache = nullptr) {

    unsigned int n = t1.length();
    unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return (m == n) ? (0 <= epsilon) : (epsilon <= INFINITY);
    }


    if(n == 1 || m == 1){
        float max_dist = -1;

        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                float d = (cache) ? cache->cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));
                if(d > max_dist){
                    max_dist = d;
                }
                if (d > epsilon){
                    return false;
                }
            }
        }

        return (epsilon <= max_dist);

    }

    // Check First and Last Pairs
    float start_distance = T::dist(t1.at(0), t2.at(0));
    float final_distance = T::dist(t1.at(n - 1), t2.at(m - 1));

    if((start_distance > epsilon) || (final_distance > epsilon)){
        return false;
    }

    // Pseudo Priority Queue
    // Cells to the right will be inserted at the front, cells to the top will be inserted at the back
    // in order to preserve a valid traversal of the free space diagram
    auto q = std::deque<std::shared_ptr<FreeSpaceCell>>();

    // Hashmap of FS cells (key =  i * m + j)
    auto fs_map = new std::shared_ptr<FreeSpaceCell>[n * m];
    Bitmask fs_visited(n * m);

    // Create FreeSpace for (0,0) and insert into queue
    FreeSpaceCell fst_cell = T::createFreeSpace(t1.at(0), t2.at(0), t1.at(1), t2.at(1), epsilon, 0, 0, cache);
    q.push_front(std::make_shared<FreeSpaceCell>(fst_cell));
    fs_map[0] = std::make_shared<FreeSpaceCell>(fst_cell);
    fs_visited.set_bit(0);

    while(!q.empty()){
        FreeSpaceCell fs = *q.front();
        q.pop_front();

        // First Cell
        if(fs.i == 0 && fs.j == 0){

            // Remove upper intervals (needed for road network free space)
            if(T::fs_intervals_per_side == 2){
                fs.setLeft(fs.getLeft().first, FreeSpaceCell::null_pair);
                fs.setBottom(fs.getBottom().first, FreeSpaceCell::null_pair);
            }

            if(!fs.in(0,LEFT)){          // Left interval does not contain (0,0)
                fs.setLeft(FreeSpaceCell::null_pair, FreeSpaceCell::null_pair);
            }

            if(!fs.in(0,BOTTOM)){          // Bottom interval does not contain (0,0)
                fs.setBottom(FreeSpaceCell::null_pair, FreeSpaceCell::null_pair);
            }

        }

        // Adjust right and top bounds of cell to be reachable
        T::propagateFreeSpace(t1.at(fs.i), t2.at(fs.j), t1.at(fs.i + 1), t2.at(fs.j + 1), fs, epsilon, cache);

        // Last Cell
        if(fs.i == (n - 2) && fs.j == (m - 2)){                                     // No path to (n-1,m-1)

            delete[] fs_map;

            if(!(fs.in(n - 1, RIGHT)) || !(fs.in(m - 1, TOP))){
                return false;
            }

            return true; // Path to (n-1,m-1) was found

        }


        // Propagate to the right
        if((fs.getRight().first != FreeSpaceCell::null_pair || fs.getRight().second != FreeSpaceCell::null_pair) && fs.j < m - 2){

            if(!fs_visited.bit_set(fs.i * m + fs.j + 1)){   // Cell to the right has not been created yet
                auto fs_r = FreeSpaceCell(fs.i, fs.j + 1, T::fs_intervals_per_side);
                fs_r.setLeft(fs.getRight().first, fs.getRight().second);

                if(fs.i == 0){
                    T::computeFreeSpaceSide(t1.at(fs.i), t2.at(fs.j + 1), t1.at(fs.i + 1), t2.at(fs.j + 2), fs_r, epsilon, BOTTOM, cache);   // If we're at the bottommost row, we need to compute the lower bound
                }

                q.push_front(std::make_shared<FreeSpaceCell>(fs_r));
                fs_map[fs.i * m + fs.j + 1] = std::make_shared<FreeSpaceCell>(fs_r);
                fs_visited.set_bit(fs.i * m + fs.j + 1);
            }else{
                fs_map[fs.i * m + fs.j + 1]->setLeft(fs.getRight().first, fs.getRight().second);
            }

        }

        // Propagate to the top
        if((fs.getTop().first != FreeSpaceCell::null_pair || fs.getTop().second != FreeSpaceCell::null_pair) && fs.i < n - 2){

            if(!fs_visited.bit_set((fs.i + 1) * m + fs.j)){   // Cell to the top has not been created yet
                auto fs_t = FreeSpaceCell(fs.i + 1, fs.j, T::fs_intervals_per_side);
                fs_t.setBottom(fs.getTop().first, fs.getTop().second);

                if(fs.j == 0){
                    T::computeFreeSpaceSide(t1.at(fs.i + 1), t2.at(fs.j), t1.at(fs.i + 2), t2.at(fs.j + 1), fs_t, epsilon, LEFT, cache);    // If we're at the leftmost column, we need to compute the left bound
                }

                q.push_back(std::make_shared<FreeSpaceCell>(fs_t));
                fs_map[(fs.i + 1) * m + fs.j] = std::make_shared<FreeSpaceCell>(fs_t);
                fs_visited.set_bit((fs.i + 1) * m + fs.j);
            }
            else{
                fs_map[(fs.i + 1) * m + fs.j]->setBottom(fs.getTop().first, fs.getTop().second);
            }
        }
    }

    delete[] fs_map;
    return false; // No path found

}

template<class T>
float frechet(Trajectory<T>& t1, Trajectory<T>& t2, std::shared_ptr<DistCache<T>> cache = NULL){

    unsigned int n = t1.length();
    unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return (m == n) ? 0 : INFINITY;
    }

    if(n == 1 || m == 1){
        float max_dist = -1;

        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                float d = (cache) ? cache->cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));
                if(d > max_dist){
                    max_dist = d;
                }
            }
        }

        return max_dist;

    }

    std::vector<float> critical_values;

    // Determine all critical values of type a (as described in the algorithm by Alt et Godau)
    // i.e. Distances between start and end points
    float start_dist = (cache) ? cache->cacheDistance(0,0) : T::dist(t1.at(0), t2.at(0));
    float end_dist = (cache) ? cache->cacheDistance(n - 1, m - 1) : T::dist(t1.at(n - 1), t2.at(m - 1));
    critical_values.push_back(start_dist);
    critical_values.push_back(end_dist);

    // Determine all critical values of type b
    // i.e. Distance for each point of t1 (t2) and each edge of t2 (t1)

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m - 1; j++){
            float d = T::dist_to_edge(t1.at(i), t2.at(j), t2.at(j + 1));
            critical_values.push_back(d);
        }
    }

    for(int j = 0; j < m; j++){
        for(int i = 0; i < n - 1; i++){
            float d = T::dist_to_edge(t2.at(j), t1.at(i), t1.at(i + 1));
            critical_values.push_back(d);
        }
    }

    // Determine all critical values of type c
    // i.e. For each point pair of t1 (t2), find the distance to the intersection of the bisector and each edge in t2 (t1)
    // (if intersection exists)

    for (int k = 0; k < n; k++){
        for(int l = k + 1; l < n; l++){
            for(int j = 0; j < m - 1; j++){
                float d = T::dist_of_bisector(t1.at(k), t1.at(l), t2.at(j), t2.at(j + 1));
                if(d >= 0){
                    critical_values.push_back(d);
                }
            }
        }
    }

    for (int k = 0; k < m; k++){
        for(int l = k + 1; l < m; l++){
            for(int i = 0; i < n - 1; i++){
                float d = T::dist_of_bisector(t2.at(k), t2.at(l), t1.at(i), t1.at(i + 1));
                if(d >= 0){
                    critical_values.push_back(d);
                }
            }
        }
    }

    // Sort and do a binary search over critical values
    // (Yes, the paper calls for a parametric search, but such a parametric search is both complex and far more inefficient in practice)
    std::sort(critical_values.begin(), critical_values.end());

    int l = 0;
    int r = critical_values.size() - 1;

    while(l < r){

        int med = l + (int)((r - l)/2);

        if(frechet_LT_Epsilon(t1,t2,critical_values[med])){
            r = med - 1;
        }else{
            l = med + 1;
        }

    }

    return critical_values[l];

}


