/*
 *      ==============================
 *          These Measures are variations only meant to be used in the Experiments
 *      ==============================
 */

#pragma once

#include "ModifiedFrechet.hpp"
#include "ModifiedHausdorff.hpp"
#include "ModifiedDTW.hpp"
#include "CombinedSimilarityMeasures.hpp"

/*
 *          ====================
 *             TIMED VARIANTS
 *          ====================
 */

template<class T>
std::tuple<double,double,double> timed_fr_hd_dtw_algorithm(Trajectory<T>& t1, Trajectory<T>& t2, float fr_eps, float hd_eps, float dtw_eps, bool use_cache = false){
    auto hd_ind = HausdorffIndexManager(t1.length(), t2.length());
    auto dtw_paths = DTWPathManager(t1.length(), t2.length());

    auto cache = (use_cache) ? std::make_shared<DistCache<T>>(std::make_shared<Trajectory<T>>(t1),std::make_shared<Trajectory<T>>(t2)) : nullptr;

    double fr_rt = 0, hd_rt = 0, dtw_rt = 0;

    auto start = std::chrono::high_resolution_clock::now();
    modified_frechet_LT_epsilon(t1,t2,fr_eps, hd_eps, hd_ind, dtw_paths, cache);
    auto end = std::chrono::high_resolution_clock::now();
    fr_rt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    if(hd_ind.decided()){
        hd_ind.lt_eps();
    }else{
        modified_hd_LT_epsilon(t1,t2,hd_eps, hd_ind, cache);
    }
    end = std::chrono::high_resolution_clock::now();
    hd_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    if(auto dtw_bound = dtw_paths.get_bound(t1.length() - 1,t2.length() - 1); dtw_bound.second != NONE){
        if(dtw_bound.second == GREATER && dtw_bound.first > dtw_eps){
            false;
        }
        else if(dtw_bound.second == LESS_EQ && dtw_bound.first <= dtw_eps){
            true;
        }
        else{
            modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
        }

    }else{
        modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
    }
    end = std::chrono::high_resolution_clock::now();
    dtw_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    return  {fr_rt, hd_rt, dtw_rt};

}


template<class T>
std::tuple<double,double,double> timed_fr_hd_dtw_similarity_join(Trajectory<T>& t1, Trajectory<T>& t2, float fr_eps, float hd_eps, float dtw_eps, COMB_MEASURE combiner = AND, bool use_cache = false){
    bool fr_LT_eps, hd_LT_eps, dtw_LT_eps;

    auto hd_ind = HausdorffIndexManager(t1.length(), t2.length());
    auto dtw_paths = DTWPathManager(t1.length(), t2.length());

    auto cache = (use_cache) ? std::make_shared<DistCache<T>>(std::make_shared<Trajectory<T>>(t1),std::make_shared<Trajectory<T>>(t2)) : nullptr;

    double fr_rt = 0, hd_rt = 0, dtw_rt = 0;

    auto start = std::chrono::high_resolution_clock::now();
    fr_LT_eps = modified_frechet_LT_epsilon(t1,t2,fr_eps, hd_eps, hd_ind, dtw_paths, cache);
    auto end = std::chrono::high_resolution_clock::now();
    fr_rt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // Early Breaks
    if(combiner == AND && !fr_LT_eps){ return {fr_rt, 0, 0};}
    else if(combiner == OR && fr_LT_eps){return {fr_rt, 0, 0};}

    start = std::chrono::high_resolution_clock::now();
    if(hd_ind.decided()){
        hd_LT_eps = hd_ind.lt_eps();
    }else{
        hd_LT_eps = modified_hd_LT_epsilon(t1,t2,hd_eps, hd_ind, cache);
    }
    end = std::chrono::high_resolution_clock::now();
    hd_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();


    // Early Breaks
    if(combiner == AND && !hd_LT_eps){return {fr_rt, hd_rt, 0};}
    else if(combiner == OR && hd_LT_eps){return {fr_rt, hd_rt, 0};}
    else if(combiner == MAJORITY && (fr_LT_eps == hd_LT_eps)){return {fr_rt, hd_rt, 0};}

    start = std::chrono::high_resolution_clock::now();
    if(auto dtw_bound = dtw_paths.get_bound(t1.length() - 1,t2.length() - 1); dtw_bound.second != NONE){
        if(dtw_bound.second == GREATER && dtw_bound.first > dtw_eps){
            dtw_LT_eps = false;
        }
        else if(dtw_bound.second == LESS_EQ && dtw_bound.first <= dtw_eps){
            dtw_LT_eps = true;
        }
        else{
            dtw_LT_eps = modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
        }

    }else{
        dtw_LT_eps = modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
    }
    end = std::chrono::high_resolution_clock::now();
    dtw_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    return {fr_rt, hd_rt, dtw_rt};

}

/*
 *          ========================
 *             FS EXPLORAL VARIANTS
 *          ========================
 */

template<class T>
int fs_exploral_fr_LT_epsilon(Trajectory<T>& t1, Trajectory<T>& t2, float fr_epsilon, float hd_epsilon, HausdorffIndexManager& hd_ind, DTWPathManager& dtw_mgr, std::shared_ptr<DistCache<T>> cache = NULL){ // Attention: Cache never truly used!

    unsigned int n = t1.length();
    unsigned int m = t2.length();

    int expl_cells = 0;

    if(n == 0 || m == 0){
        return 0;
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
                    return 0;
                }
            }
        }

        return 0;

    }

    // Check First and Last Pairs
    float start_distance = T::dist(t1.at(0), t2.at(0));
    float final_distance = T::dist(t1.at(n - 1), t2.at(m - 1));

    if((start_distance > fr_epsilon) || (final_distance > fr_epsilon)){
        return 0;
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
    FreeSpaceCell fst_cell = T::createFreeSpace(t1.at(0), t2.at(0), t1.at(1), t2.at(1), fr_epsilon, 0, 0);
    q.push_front(std::make_shared<FreeSpaceCell>(fst_cell));
    fs_map[0] = std::make_shared<FreeSpaceCell>(fst_cell);
    fs_visited.set_bit(0);

    while(!q.empty()){
        FreeSpaceCell fs = *q.front();
        q.pop_front();

        expl_cells++;

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
        T::propagateFreeSpace(t1.at(fs.i), t2.at(fs.j), t1.at(fs.i + 1), t2.at(fs.j + 1), fs, fr_epsilon);

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

                return expl_cells;
            }

            return expl_cells;                                                                  // Path to (n-1,m-1) was found

        }


        // Propagate to the right
        if((fs.getRight().first != FreeSpaceCell::null_pair || fs.getRight().second != FreeSpaceCell::null_pair) && fs.j < m - 2){

            if(!fs_visited.bit_set(fs.i * m + fs.j + 1)){                           // Cell to the right has not been created yet

                auto fs_r = FreeSpaceCell(fs.i, fs.j + 1, T::fs_intervals_per_side);
                fs_r.setLeft(fs.getRight().first, fs.getRight().second);

                if(fs.i == 0){
                    T::computeFreeSpaceSide(t1.at(fs.i), t2.at(fs.j + 1), t1.at(fs.i + 1), t2.at(fs.j + 2), fs_r, fr_epsilon, BOTTOM);   // If we're at the bottommost row, we need to compute the lower bound
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
                    T::computeFreeSpaceSide(t1.at(fs.i + 1), t2.at(fs.j), t1.at(fs.i + 2), t2.at(fs.j + 1), fs_t, fr_epsilon, LEFT);    // If we're at the leftmost column, we need to compute the left bound
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
    return expl_cells; // No path found

}

template<class T>
std::tuple<int,double,double, double> timed_exploral_fr_hd_dtw_algorithm(Trajectory<T>& t1, Trajectory<T>& t2, float fr_eps, float hd_eps, float dtw_eps, bool use_cache = false){
    auto hd_ind = HausdorffIndexManager(t1.length(), t2.length());
    auto dtw_paths = DTWPathManager(t1.length(), t2.length());

    auto cache = (use_cache) ? std::make_shared<DistCache<T>>(std::make_shared<Trajectory<T>>(t1),std::make_shared<Trajectory<T>>(t2)) : nullptr;

    double fr_rt = 0, hd_rt = 0, dtw_rt = 0;

    auto start = std::chrono::high_resolution_clock::now();
    int expl_fs = fs_exploral_fr_LT_epsilon(t1,t2,fr_eps, hd_eps, hd_ind, dtw_paths, cache);
    auto end = std::chrono::high_resolution_clock::now();
    fr_rt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    if(hd_ind.decided()){
        hd_ind.lt_eps();
    }else{
        modified_hd_LT_epsilon(t1,t2,hd_eps, hd_ind, cache);
    }
    end = std::chrono::high_resolution_clock::now();
    hd_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    if(auto dtw_bound = dtw_paths.get_bound(t1.length() - 1,t2.length() - 1); dtw_bound.second != NONE){
        if(dtw_bound.second == GREATER && dtw_bound.first > dtw_eps){
            false;
        }
        else if(dtw_bound.second == LESS_EQ && dtw_bound.first <= dtw_eps){
            true;
        }
        else{
            modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
        }

    }else{
        modified_dtw_LT_epsilon(t1,t2,dtw_eps,cache);
    }
    end = std::chrono::high_resolution_clock::now();
    dtw_rt = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    return  {expl_fs, fr_rt, hd_rt, dtw_rt};

}