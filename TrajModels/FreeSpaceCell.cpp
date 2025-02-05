#include "FreeSpaceCell.hpp"

std::pair<fpair, fpair> FreeSpaceCell::getLeft() {
    return {left_low, left_up};
}

std::pair<fpair, fpair> FreeSpaceCell::getBottom() {
    return {bottom_low, bottom_up};
}

std::pair<fpair, fpair> FreeSpaceCell::getRight() {
    return {right_low, right_up};
}

std::pair<fpair, fpair> FreeSpaceCell::getTop() {
    return {top_low, top_up};
}

std::string FreeSpaceCell::to_string() {

    return "Left\t" + interval_to_string(LEFT) + "\nBottom\t" + interval_to_string(BOTTOM) + "\nRight\t" +
           interval_to_string(RIGHT) + "\nTop\t" + interval_to_string(TOP) + "\n";
}

std::string FreeSpaceCell::interval_to_string(side s) {

    fpair interval_low;
    fpair interval_up;

    switch(s){
        case LEFT:
            interval_low = left_low;
            interval_up = left_up;
            break;
        case RIGHT:
            interval_low = right_low;
            interval_up = right_up;
            break;
        case TOP:
            interval_low = top_low;
            interval_up = top_up;
            break;
        case BOTTOM:
            interval_low = bottom_low;
            interval_up = bottom_up;
            break;
    }

    if(interval_low == null_pair && interval_up == null_pair){
        return "<empty>";
    }
    else {
        std::string ret;
        ret += (interval_low == null_pair) ? "<empty> | " : (std::to_string(interval_low.first) + " " + std::to_string(interval_low.second));
        ret += (interval_up == null_pair) ? "   <empty>" : "   " + (std::to_string(interval_up.first) + " " + std::to_string(interval_up.second));
        return ret;
    }



}


/*
        Checks if a floating point value is contained on a side of the free space diagram
 */
bool FreeSpaceCell::in(float x, side s) {

    fpair interval_low;
    fpair interval_up;

    switch(s){
        case LEFT:
            interval_low = left_low;
            interval_up = left_up;
            break;
        case RIGHT:
            interval_low = right_low;
            interval_up = right_up;
            break;
        case TOP:
            interval_low = top_low;
            interval_up = top_up;
            break;
        case BOTTOM:
            interval_low = bottom_low;
            interval_up = bottom_up;
            break;
    }

    bool inside = false;

    if(interval_low.first != -1 || interval_low.second != -1){
        inside = (interval_low.first <= x && interval_low.second >= x);
    }
    if(interval_up.first != -1 || interval_up.second != -1){
        inside = inside || (interval_up.first <= x && interval_up.second >= x);
    }

    return inside;
}

void FreeSpaceCell::setLeft(fpair left_low) {
    this->left_low = left_low;
}

void FreeSpaceCell::setLeft(fpair left_low, fpair left_up) {
    this->left_low = left_low;
    this->left_up = left_up;
}

void FreeSpaceCell::setRight(fpair right_low) {
    this->right_low = right_low;
}

void FreeSpaceCell::setRight(fpair right_low, fpair right_up) {
    this->right_low = right_low;
    this->right_up  = right_up;
}

void FreeSpaceCell::setBottom(fpair bottom_low) {
    this->bottom_low = bottom_low;
}

void FreeSpaceCell::setBottom(fpair bottom_low, fpair bottom_up) {
    this->bottom_low = bottom_low;
    this->bottom_up  = bottom_up;
}

void FreeSpaceCell::setTop(fpair top_low) {
    this->top_low = top_low;
}

void FreeSpaceCell::setTop(fpair top_low, fpair top_up) {
    this->top_low = top_low;
    this->top_up  = top_up;
}

bool FreeSpaceCell::bl_reachable() {
    return in(i,LEFT) && in(j,BOTTOM);
}

bool FreeSpaceCell::tl_reachable() {
    return in(i + 1,LEFT) && in(j,TOP);
}

bool FreeSpaceCell::tr_reachable() {
    return in(i + 1,RIGHT) && in(j + 1,TOP);
}

bool FreeSpaceCell::br_reachable() {
    return in(i,RIGHT) && in(j + 1,BOTTOM);
}


