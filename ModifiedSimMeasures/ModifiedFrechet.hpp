#pragma once
#include <deque>
#include <algorithm>
#include <unordered_map>
#include <cmath>

#include "../TrajModels/Trajectory.hpp"
#include "../TrajModels/FreeSpaceCell.hpp"
#include "../DistCache/DistCache.hpp"

#include "HausdorffIndexManager.hpp"
#include "DTWPathManager.hpp"

template<class T>
void checkPoint(int i, int j, bool reachable, int n, int m, float fr_epsilon, float hd_epsilon, HausdorffIndexManager& hd_ind, DTWPathManager& dtw_mgr, int* rows_grey, int* cols_grey, Bitmask& rows_seen, Bitmask& columns_seen){

    dtw_mgr.store_bound(i,j,reachable,fr_epsilon);

    if(reachable){
        if(fr_epsilon > hd_epsilon){
            hd_ind.insert_unprioritized(i,j);
        }else{
            hd_ind.insert_blacklisted(i,j);
        }

    }
    else{

        bool curr_row_seen = rows_seen.bit_set(i);
        bool curr_col_seen = columns_seen.bit_set(j);

        if(curr_row_seen && !curr_col_seen){
            rows_grey[i]++;

            if(rows_grey[i] == m){ // entire row grey
                (fr_epsilon >= hd_epsilon) ? hd_ind.decide_hausdorff(false) : hd_ind.insert_prioritized(i,true);
            }

        }else if (!curr_row_seen){
            rows_grey[i] = 1;
        }

        if(curr_col_seen && !curr_row_seen){
            cols_grey[j]++;

            if(cols_grey[j] == n){ // entire column grey
                (fr_epsilon >= hd_epsilon) ? hd_ind.decide_hausdorff(false) : hd_ind.insert_prioritized(j,false);
            }
        }
        else if(!curr_col_seen){
            cols_grey[j] = 1;
        }

        if(!curr_row_seen) {rows_seen.set_bit(i);}
        if(!curr_col_seen) {columns_seen.set_bit(j);}

    }

}

template<class T>
bool modified_frechet_LT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float fr_epsilon, float hd_epsilon, HausdorffIndexManager& hd_ind, DTWPathManager& dtw_mgr, std::shared_ptr<DistCache<T>> cache = NULL){ // Attention: Cache never truly used!

    unsigned int n = t1.length();
    unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return (m == n) ? (0 <= fr_epsilon) : (fr_epsilon <= INFINITY);
    }

    if(n == 1 || m == 1){
        float max_dist = -1;

        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                float d = (cache) ? cache->cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));
                if(d > max_dist){
                    max_dist = d;
                }
                if (d > fr_epsilon){
                    return false;
                }
            }
        }

        return (fr_epsilon <= max_dist);

    }

    // Check First and Last Pairs
    float start_distance = T::dist(t1.at(0), t2.at(0));
    float final_distance = T::dist(t1.at(n - 1), t2.at(m - 1));

    if((start_distance > fr_epsilon) || (final_distance > fr_epsilon)){
        return false;
    }

    // Pseudo Priority Queue
    // Cells to the right will be inserted at the front, cells to the top will be inserted at the back
    // in order to preserve a valid traversal of the free space diagram
    auto q = std::deque<std::shared_ptr<FreeSpaceCell>>();

    // Hashmap of FS cells (unordered map) (key =  i * m + j)

    auto fs_map = new std::shared_ptr<FreeSpaceCell>[n * m];
    Bitmask fs_visited(n * m);

    // Hashmaps keeping track of how many grey points have been found in each column/row
    int* rows_grey = new int[n];
    int* columns_grey = new int[m];
    auto rows_seen = Bitmask(n);
    auto columns_seen = Bitmask(m);

    // Create FreeSpace for (0,0) and insert into queue
    FreeSpaceCell fst_cell = T::createFreeSpace(t1.at(0), t2.at(0), t1.at(1), t2.at(1), fr_epsilon, 0, 0, cache);
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

            // check (0,0)
            checkPoint<T>(0,0,fs.bl_reachable(), n, m, fr_epsilon, hd_epsilon, hd_ind, dtw_mgr, rows_grey, columns_grey, rows_seen, columns_seen);

        }

        // Adjust right and top bounds of cell to be reachable
        T::propagateFreeSpace(t1.at(fs.i), t2.at(fs.j), t1.at(fs.i + 1), t2.at(fs.j + 1), fs, fr_epsilon, cache);

        // Check leftmost and bottommost points (for first row/column)
        if(fs.j == 0){
            checkPoint<T>(fs.i + 1, fs.j, fs.tl_reachable(), n, m, fr_epsilon, hd_epsilon, hd_ind, dtw_mgr, rows_grey, columns_grey, rows_seen, columns_seen);
        }
        else if(fs.i == 0){
            checkPoint<T>(fs.i, fs.j + 1, fs.br_reachable(), n, m, fr_epsilon, hd_epsilon, hd_ind, dtw_mgr, rows_grey, columns_grey, rows_seen, columns_seen);
        }

        // If left/diagonal/bottom cell has never been visited, check the top left/bottom left/bottom right points-8´ß
        // These cases are distinct from the above ones, because the indexing of our fs_map runs into problems for j < 0 or i < 0

        if(fs.j != 0 && (!fs_visited.bit_set(fs.i * m + fs.j - 1))){
            checkPoint<T>(fs.i + 1, fs.j, fs.tl_reachable(), n, m, fr_epsilon, hd_epsilon, hd_ind, dtw_mgr, rows_grey, columns_grey, rows_seen, columns_seen);

            fs_visited.set_bit(fs.i * m + fs.j - 1);     // To avoid checking the same point multiple times, we mark this FS cell as seen
                                                            // This will not affect the course of the normal algorithm since we never look back at left/bottom/bottom left cells while traversing Free Space
        }
        if(fs.i != 0 && fs.j != 0 && (!fs_visited.bit_set((fs.i - 1) * m + fs.j - 1))){
            checkPoint<T>(fs.i, fs.j, fs.bl_reachable(), n, m, fr_epsilon, hd_epsilon, hd_ind, dtw_mgr, rows_grey, columns_grey, rows_seen, columns_seen);
            fs_visited.set_bit((fs.i - 1) * m + fs.j - 1);
        }
        if(fs.i != 0 && (!fs_visited.bit_set((fs.i - 1) * m + fs.j))){
            checkPoint<T>(fs.i,fs.j + 1, fs.br_reachable(), n, m, fr_epsilon, hd_epsilon, hd_ind, dtw_mgr, rows_grey, columns_grey, rows_seen, columns_seen);
            fs_visited.set_bit((fs.i - 1) * m + fs.j);
        }

        // Check top right point
        checkPoint<T>(fs.i + 1, fs.j + 1, fs.tr_reachable(), n, m, fr_epsilon, hd_epsilon, hd_ind, dtw_mgr, rows_grey, columns_grey, rows_seen, columns_seen);

        // Last Cell
        if(fs.i == (n - 2) && fs.j == (m - 2)){

            delete[] fs_map;
            delete[] rows_grey;
            delete[] columns_grey;

            if(!(fs.in(n - 1, RIGHT)) || !(fs.in(m - 1, TOP))){             // No path to (n-1,m-1)

                    return false;
            }

            return true;                                                                    // Path to (n-1,m-1) was found

        }


        // Propagate to the right
        if((fs.getRight().first != FreeSpaceCell::null_pair || fs.getRight().second != FreeSpaceCell::null_pair) && fs.j < m - 2){

            if(!fs_visited.bit_set(fs.i * m + fs.j + 1)){                           // Cell to the right has not been created yet

                auto fs_r = FreeSpaceCell(fs.i, fs.j + 1, T::fs_intervals_per_side);
                fs_r.setLeft(fs.getRight().first, fs.getRight().second);

                if(fs.i == 0){
                    T::computeFreeSpaceSide(t1.at(fs.i), t2.at(fs.j + 1), t1.at(fs.i + 1), t2.at(fs.j + 2), fs_r, fr_epsilon, BOTTOM, cache);   // If we're at the bottommost row, we need to compute the lower bound
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

            if(!fs_visited.bit_set((fs.i + 1) * m + fs.j)){                             // Cell to the top has not been created yet

                auto fs_t = FreeSpaceCell(fs.i + 1, fs.j, T::fs_intervals_per_side);
                fs_t.setBottom(fs.getTop().first, fs.getTop().second);

                if(fs.j == 0){
                    T::computeFreeSpaceSide(t1.at(fs.i + 1), t2.at(fs.j), t1.at(fs.i + 2), t2.at(fs.j + 1), fs_t, fr_epsilon, LEFT, cache);    // If we're at the leftmost column, we need to compute the left bound
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
    delete[] rows_grey;
    delete[] columns_grey;
    return false; // No path found


}


template<class T>
bool modified_discrete_frechet_LT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float fr_epsilon, float hd_epsilon, HausdorffIndexManager& hd_ind, DTWPathManager& dtw_mgr, std::shared_ptr<DistCache<T>> cache = NULL){
    const unsigned int n = t1.length();
    const unsigned int m = t2.length();

    if(n == 0 || m == 0){
        return (m == n);
    }

    auto *dFr = new float[n*m];

    bool curr_diag = true;				// True if all values are > epsilon
    bool prev_diag = false;

    int* rows_grey = new int[n];
    int* columns_grey = new int[m];
    auto rows_seen = Bitmask(n);
    auto columns_seen = Bitmask(m);

    // Traverse Matrix in diagonal fashion
    for (int k = 0; k < (n + m - 1); k++){

        int i;
        int j;

        // Early Break
        if(prev_diag && curr_diag){
            delete[] dFr;
            delete[] rows_grey;
            delete[] columns_grey;
            return false;
        }

        if(k < n){   // for diagonals before (and including) main diagonal
            i = k;
            j = 0;
        }
        else{
            i = n - 1;
            j = k - (n - 1);
        }

        prev_diag = curr_diag;
        curr_diag = true;

        // Traverse diagonal
        while(i >= 0 && j < m){

            float d = cache ? (*cache).cacheDistance(i,j) : T::dist(t1.at(i), t2.at(j));

            checkPoint<T>(i,j, (d <= fr_epsilon), n, m, fr_epsilon, hd_epsilon, hd_ind, dtw_mgr, rows_grey, columns_grey, rows_seen, columns_seen);

            if (i == 0 && j == 0) {
                dFr[i * m + j] = d;
            }
            else if (i == 0) {
                dFr[i * m + j] = std::max(d, dFr[i * m + (j - 1)]);
            }
            else if (j == 0) {
                dFr[i * m + j] = std::max(d, dFr[(i - 1) * m + j]);
            }
            else {
                dFr[i * m + j] = std::max(d, std::min({ dFr[i * m + (j - 1)], dFr[(i - 1) * m + j], dFr[(i - 1) * m + (j - 1)] }));
            }

            if(dFr[i * m + j] < fr_epsilon){
                curr_diag = false;
            }

            i--;
            j++;
        }

    }

    float disc_fr = dFr[(n - 1) * m + (m - 1)];
    delete[] dFr;
    delete[] rows_grey;
    delete[] columns_grey;
    return (disc_fr < fr_epsilon);
}