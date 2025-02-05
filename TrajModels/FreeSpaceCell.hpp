#pragma once 
#include <memory>
#include <utility>
#include <string>

enum side  {LEFT, TOP, RIGHT, BOTTOM};

typedef std::pair<float, float> fpair;

/*
    Models the Free Space cells used by the Alt and Godau algorithm for the Frechet distance
 */
class FreeSpaceCell {

    private: 
        int intervals_per_side;

        fpair left_low = null_pair;
        fpair left_up = null_pair;
        fpair bottom_low = null_pair;
        fpair bottom_up = null_pair;
        fpair top_low = null_pair;
        fpair top_up = null_pair;
        fpair right_low = null_pair;
        fpair right_up = null_pair;


public:

    inline static fpair null_pair = std::pair<float,float>(-1,-1);

    int i;
        int j;

        FreeSpaceCell(int i, int j, int intervals_per_side) : i(i), j(j), intervals_per_side(intervals_per_side) {};
        FreeSpaceCell(int i, int j, int intervals_per_side, fpair left_low, fpair bottom_low, fpair right_low, fpair top_low, fpair left_up = null_pair, fpair bottom_up = null_pair, fpair right_up = null_pair, fpair top_up = null_pair)
                      : i(i), j(j), intervals_per_side(intervals_per_side), left_low(std::move(left_low)), bottom_low(std::move(bottom_low)), right_low(std::move(right_low)), top_low(std::move(top_low)), left_up(std::move(left_up)), bottom_up(std::move(bottom_up)), right_up(std::move(right_up)), top_up(std::move(top_up)){};

        std::pair<fpair, fpair> getLeft();
        std::pair<fpair, fpair> getBottom();
        std::pair<fpair, fpair> getRight();
        std::pair<fpair, fpair> getTop();

        void setLeft(fpair left_low);
        void setLeft(fpair left_low, fpair left_up);

        void setRight(fpair right_low);
        void setRight(fpair right_low, fpair right_up);

        void setBottom(fpair bottom_low);
        void setBottom(fpair bottom_low, fpair bottom_up);

        void setTop(fpair top_low);
        void setTop(fpair top_low, fpair top_up);

        std::string to_string();
        std::string interval_to_string(side s);

        bool in(float x, side s);

        // Corner points
        bool bl_reachable();
        bool tl_reachable();
        bool tr_reachable();
        bool br_reachable();

};