#pragma once
#include <memory>
#include <cmath>
#include "../Utils/Bitmask.hpp"

enum BOUND {LESS_EQ, GREATER, NONE};

class DTWPathManager{
    private:

        float * path_val;
        BOUND * path_bound;
        Bitmask cached;

        int col_nr = -1;

        //std::unordered_map<int,int> cached_diags;

        bool contains(int i, int j){
            return cached.bit_set(i * col_nr + j);
        }

        void insert(int i, int j, BOUND b, float value){
            cached.set_bit(i * col_nr + j);
            path_bound[i * col_nr + j] = b;
            path_val[i * col_nr + j] = value;
        }

        static float lowest_bound(BOUND b1, float v1, BOUND b2, float v2, BOUND b3, float v3){
            float min = INFINITY;

            // First look for all LESS_EQ bounds
            if(b1 == LESS_EQ) {min = v1;}
            if(b2 == LESS_EQ && v2 < min) {min = v2;}
            if(b3 == LESS_EQ && v3 < min) {min = v3;}
            if(min != INFINITY) {return min;}

            // If there are no LESS_EQ bounds, then look for all GREATER bounds
            if(b1 == GREATER) {min = v1;}
            if(b2 == GREATER && v2 < min) {min = v2;}
            if(b3 == GREATER && v3 < min) {min = v3;}
            return min;

        }

    public:

        DTWPathManager(int n, int m) : col_nr(m), cached(Bitmask(n * m)){
            path_bound = new BOUND[n * m];
            path_val = new float[n * m];
        }

        ~DTWPathManager(){
            delete[] path_val;
            delete[] path_bound;
        }

        void store_bound(int i, int j, bool reachable, float epsilon){

            if(i == 0 && j == 0){
                (reachable) ? insert(0,0,LESS_EQ, epsilon) : insert(0,0,GREATER,epsilon);
                return;
            }

            else if(i == 0){
                auto left_bound = (contains(i,j-1)) ? path_bound[i * col_nr + j - 1] : NONE;
                auto left_val = (contains(i,j-1)) ? path_val[i * col_nr + j - 1] : INFINITY;

                if(left_bound == GREATER){
                    (reachable) ? insert(i,j,GREATER,left_val) : insert(i,j,GREATER, left_val + epsilon);
                    return;
                }
                else if(left_bound == LESS_EQ){
                    (reachable) ? insert(i,j,LESS_EQ, left_val + epsilon) : insert(i,j,GREATER,epsilon);
                    return;
                }
                else{
                    (reachable) ? insert(i,j,NONE,INFINITY) : insert(i,j,GREATER,epsilon);
                    return;
                }
            }

            else if(j == 0){

                auto bottom_bound = (contains(i - 1, j)) ? path_bound[(i - 1) * col_nr + j] : NONE;
                auto bottom_val = (contains(i - 1, j)) ? path_val[(i - 1) * col_nr + j] : INFINITY;

                if(bottom_bound == GREATER){
                    (reachable) ? insert(i,j,GREATER,bottom_val) : insert(i,j,GREATER, bottom_val + epsilon);
                    return;
                }
                else if(bottom_bound == LESS_EQ){
                    (reachable) ? insert(i,j,LESS_EQ, bottom_val + epsilon) : insert(i,j,GREATER,epsilon);
                    return;
                }
                else{
                    (reachable) ? insert(i,j,NONE,INFINITY) : insert(i,j,GREATER,epsilon);
                    return;
                }

            }

            // For the reasonings of the bounds that follow, see the thesis


            BOUND left_bound = NONE, bottom_bound = NONE, diag_bound = NONE;
            float left_val = INFINITY, bottom_val = INFINITY, diag_val = INFINITY;

            if(contains(i,j - 1)) {
                int index = i * col_nr + j - 1;
                left_bound = path_bound[index];
                left_val = path_val[index];
            }
            if(contains(i - 1, j)){
                int index = (i - 1) * col_nr + j;
                bottom_bound = path_bound[index];
                bottom_val = path_val[index];
            }
            if(contains(i - 1, j - 1)){
                int index = (i - 1) * col_nr + j - 1;
                diag_bound = path_bound[index];
                diag_val = path_val[index];
            }

            if(left_bound == LESS_EQ || bottom_bound == LESS_EQ || diag_bound == LESS_EQ){

                if(!reachable){
                    insert(i,j,GREATER,epsilon);
                    return;
                }

                float k = lowest_bound(left_bound, left_val, bottom_bound, bottom_val, diag_bound, diag_val);
                insert(i,j,LESS_EQ, k + epsilon);
                return;

            }
            else if(left_bound == NONE || bottom_bound == NONE || diag_bound == NONE){

                (reachable) ? insert(i,j,NONE,INFINITY) : insert(i,j,GREATER,epsilon);
                return;

            }
            else{
                float k = lowest_bound(left_bound, left_val, bottom_bound, bottom_val, diag_bound, diag_val);
                (reachable) ? insert(i,j,GREATER, k) : insert(i,j,GREATER, k + epsilon);
                return;
            }

        }

        std::pair<float,BOUND> get_bound(int i, int j){
            if(!contains(i,j)){
                return {INFINITY, NONE};
            }

            return {path_val[i * col_nr + j],path_bound[i * col_nr + j]};
        }

};