#include <iostream>
#include "RoadNetworkPoint.hpp"

std::shared_ptr<Graph> RoadNetworkPoint::g = nullptr;
std::shared_ptr<PrunedHighwayLabeling> RoadNetworkPoint::phl = nullptr;


/*
    Assumes that the Graph has a weightmap property already and that i and j stem from the same graph
    Distance is queries using PHL labels
 */
float RoadNetworkPoint::dist(RoadNetworkPoint &i, RoadNetworkPoint &j) {

    return phl->Query(i.node_id, j.node_id);

}

std::string RoadNetworkPoint::to_string() {
    return "Node " + std::to_string(node_id);
}

/*
    Creates the Free Space cell [p1,p2] x [q1,q2] (following the approach of https://doi.org/10.1007/978-3-642-24983-9_7)
 */
FreeSpaceCell RoadNetworkPoint::createFreeSpace(RoadNetworkPoint &p1, RoadNetworkPoint &q1, RoadNetworkPoint &p2, RoadNetworkPoint &q2, float epsilon, int i, int j, std::shared_ptr<DistCache<RoadNetworkPoint>> cache) {

    // Node distances
    float dist_p1_q1 = (cache) ? cache->cacheDistance(i, j) : RoadNetworkPoint::dist(p1,q1);
    float dist_p1_q2 = (cache) ? cache->cacheDistance(i, j + 1) : RoadNetworkPoint::dist(p1,q2);
    float dist_p2_q1 = (cache) ? cache->cacheDistance(i + 1, j) : RoadNetworkPoint::dist(p2,q1);
    float dist_p2_q2 = (cache) ? cache->cacheDistance(i + 1, j + 1) : RoadNetworkPoint::dist(p2,q2);

    // Edge costs
    typedef boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
    edge_descriptor p = boost::edge(p1.node_id, p2.node_id, *g).first;
    edge_descriptor q = boost::edge(q1.node_id, q2.node_id, *g).first;

    float cost_p = boost::get(boost::edge_weight, *g, p);
    float cost_q = boost::get(boost::edge_weight, *g, q);

    /*

     We need to calculate intersections for 4 lines. To do this we need to get from the analytical equation to the form
     y = mx + b

     As an example: L1 has the analytical form
     dist(p1,q1) + length(p_t, p_1) + length(q_1, q_s) = epsilon
     ===
     dist(p1,q1) + cost(p1,p2) * (s - i) + cost(q1,q2) * (t - j) = epsilon
     ===
     s = - (dist(p1,q1) - i * cost(p1,p2) + (-j) * cost(q1,q2) - epsilon + cost(q1,q2) * t) / cost(p1,p2)
     */

    float l1_slope = -cost_q / cost_p;
    float l1_intercept = -(dist_p1_q1 - cost_p * i - cost_q * j - epsilon)/cost_p;

    float l2_slope = cost_q / cost_p;
    float l2_intercept = (dist_p2_q1 + cost_p * (i + 1) - cost_q * j - epsilon)/cost_p;

    float l3_slope = -cost_q / cost_p;
    float l3_intercept = (dist_p2_q2 + cost_p * (i + 1) + cost_q * (j + 1) - epsilon)/cost_p;

    float l4_slope = cost_q / cost_p;
    float l4_intercept = -(dist_p1_q2 - cost_p * i + cost_q * (j + 1) - epsilon)/cost_p;

    // Left Boundary
    auto left = line_line_intersection(l1_slope, l1_intercept, l2_slope, l2_intercept, i, j, true, false);

    // Right Boundary
    auto right = line_line_intersection(l3_slope, l3_intercept, l4_slope, l4_intercept, i, j, true, true);

    // Top Boundary
    auto top = line_line_intersection(l2_slope, l2_intercept, l3_slope, l3_intercept, i, j, false, true);

    // Bottom Boundary
    auto bottom = line_line_intersection(l1_slope, l1_intercept, l4_slope, l4_intercept, i, j, false, false);

    return {i, j, fs_intervals_per_side, left.first, bottom.first, right.first, top.first, left.second, bottom.second, right.second, top.second};

}

std::pair<fpair,fpair> RoadNetworkPoint::line_line_intersection(float slope1, float intercept1, float slope2, float intercept2, int i, int j, bool vertical, bool upper) {

    float intersect_x = (intercept2 - intercept1) / (slope1 - slope2);
    float intersect_y = slope1 * intersect_x + intercept1;

    if(vertical){ // left/right boundaries

        if( (!upper && intersect_x <= j) || (upper && intersect_x >= j + 1)){ // intersection point to the left/right of boundary --> two intervals

            float y_1 = (upper) ? slope1 * (j + 1) + intercept1 : slope1 * j + intercept1;
            float y_2 = (upper) ? slope2 * (j + 1) + intercept2 : slope2 * j + intercept2;

            auto low = (std::min(y_1,y_2) <= i) ? FreeSpaceCell::null_pair : std::pair<float,float>(i, std::min(y_1, y_2));
            auto high = (std::max(y_1,y_2) >= i + 1) ? FreeSpaceCell::null_pair : std::pair<float,float>(std::max(y_1, y_2), i + 1);

            return {low, high};

        }
        else{  // intersection point to the right/left of boundary --> one interval

            auto interval = std::pair<float,float>(i, i + 1);
            return {interval, FreeSpaceCell::null_pair};

        }
    }
    else{ // bottom/top boundaries

        if( (!upper && intersect_y <= i) || (upper && intersect_y >= i + 1)){ // intersection point under/over boundary --> two intervals

            float x_1 = (upper) ? (i + 1 - intercept1) / slope1 : (i - intercept1) / slope1;
            float x_2 = (upper) ? (i + 1 - intercept2) / slope2 : (i - intercept2) / slope2;

            auto low = (std::min(x_1,x_2) <= j) ? FreeSpaceCell::null_pair : std::pair<float,float>(j, std::min(x_1,x_2));
            auto high = (std::max(x_1,x_2) >= j + 1) ? FreeSpaceCell::null_pair : std::pair<float,float>(std::max(x_1,x_2), j + 1);

            return {low, high};

        }
        else{ // intersection point to the top/bottom of boundary --> one interval

            auto interval = std::pair<float,float>(j, j + 1);
            return {fpair(j, j + 1),FreeSpaceCell::null_pair };
        }
    }


}

/*
    Adjusts the right and top boundaries of a given Free Space cell to refer to the reachable (and not just the free) space
 */
void RoadNetworkPoint::propagateFreeSpace(RoadNetworkPoint &p1, RoadNetworkPoint &q1, RoadNetworkPoint &p2, RoadNetworkPoint &q2, FreeSpaceCell &cell, float epsilon, std::shared_ptr<DistCache<RoadNetworkPoint>> cache) {

    // Create right and top boundaries (if they do not exist already)
    if(cell.getRight().first == FreeSpaceCell::null_pair || cell.getRight().second == FreeSpaceCell::null_pair){

        // Get costs/distances
        typedef boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
        edge_descriptor p = boost::edge(p1.node_id, p2.node_id, *g).first;
        edge_descriptor q = boost::edge(q1.node_id, q2.node_id, *g).first;

        float dist_p1_q2 = (cache) ? cache->cacheDistance(cell.i, cell.j + 1) : RoadNetworkPoint::dist(p1,q2);
        float dist_p2_q2 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j + 1) : RoadNetworkPoint::dist(p2,q2);

        float cost_p = boost::get(boost::edge_weight, *g, p);
        float cost_q = boost::get(boost::edge_weight, *g, q);

        // Define Lines
        float l3_slope = -cost_q / cost_p;
        float l3_intercept = (dist_p2_q2 + cost_p * (cell.i + 1) + cost_q * (cell.j + 1) - epsilon)/cost_p;

        float l4_slope = cost_q / cost_p;
        float l4_intercept = -(dist_p1_q2 - cost_p * cell.i + cost_q * (cell.j + 1) - epsilon)/cost_p;

        auto intervals = line_line_intersection(l3_slope, l3_intercept, l4_slope, l4_intercept, cell.i, cell.j, true, true);

        cell.setRight(intervals.first, intervals.second);
    }

    if(cell.getTop().first == FreeSpaceCell::null_pair || cell.getTop().second == FreeSpaceCell::null_pair){

        // Get costs/distances
        typedef boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
        edge_descriptor p = boost::edge(p1.node_id, p2.node_id, *g).first;
        edge_descriptor q = boost::edge(q1.node_id, q2.node_id, *g).first;

        float dist_p2_q1 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j) : RoadNetworkPoint::dist(p2,q1);
        float dist_p2_q2 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j + 1) : RoadNetworkPoint::dist(p2,q2);
        float cost_p = boost::get(boost::edge_weight, *g, p);
        float cost_q = boost::get(boost::edge_weight, *g, q);

        // Define Lines
        float l2_slope = cost_q / cost_p;
        float l2_intercept = (dist_p2_q1 + cost_p * (cell.i + 1) - cost_q * cell.j - epsilon)/cost_p;

        float l3_slope = -cost_q / cost_p;
        float l3_intercept = (dist_p2_q2 + cost_p * (cell.i + 1) + cost_q * (cell.j + 1) - epsilon)/cost_p;

        auto intervals = line_line_intersection(l2_slope, l2_intercept, l3_slope, l3_intercept, cell.i, cell.j, false, true);
        cell.setTop(intervals.first, intervals.second);
    }

    auto old_top= cell.getTop();  // needed for later check

    // For each interval of the top and right, we investigate its reachability

    // If lower top interval exist, its always reachable (no computation necessary)

    // Upper top interval
    if(cell.getTop().second != FreeSpaceCell::null_pair){

        // if right intersection point is to the right of border, then unreachable
        if(cell.getRight().second == FreeSpaceCell::null_pair && cell.getRight().first != std::pair<float,float>(cell.i, cell.i + 1)){
            cell.setTop(cell.getTop().first, FreeSpaceCell::null_pair);
        }
        else{ // reachable from [right intersection, right border]

            // compute right intersection
            // Get costs/distances
            typedef boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
            edge_descriptor p = boost::edge(p1.node_id, p2.node_id, *g).first;
            edge_descriptor q = boost::edge(q1.node_id, q2.node_id, *g).first;

            float dist_p1_q2 = (cache) ? cache->cacheDistance(cell.i, cell.j + 1) : RoadNetworkPoint::dist(p1,q2);
            float dist_p2_q2 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j + 1) : RoadNetworkPoint::dist(p2,q2);

            float cost_p = boost::get(boost::edge_weight, *g, p);
            float cost_q = boost::get(boost::edge_weight, *g, q);

            // Define Lines

            float l3_slope = -cost_q / cost_p;
            float l3_intercept = (dist_p2_q2 + cost_p * (cell.i + 1) + cost_q * (cell.j + 1) - epsilon)/cost_p;

            float l4_slope = cost_q / cost_p;
            float l4_intercept = -(dist_p1_q2 - cost_p * cell.i + cost_q * (cell.j + 1) - epsilon)/cost_p;

            float intersection_x = (l4_intercept - l3_intercept) / (l3_slope - l4_slope);

            auto new_upper_top = (intersection_x >= cell.j + 1) ? FreeSpaceCell::null_pair : std::pair<float,float>(intersection_x, cell.j + 1);
            cell.setTop(cell.getTop().first, new_upper_top);
        }

    }

    // If lower right interval exists, its always reachable (no computation necessary)

    // Upper right interval
    if(cell.getRight().second != FreeSpaceCell::null_pair){

        // if upper intersection point to the top of top border, then unreachable
        if(old_top.second == FreeSpaceCell::null_pair && old_top.first != std::pair<float,float>(cell.j,cell.j + 1)){
            cell.setRight(cell.getRight().first, FreeSpaceCell::null_pair);
        }
        else{ // reachable from [upper intersection, upper border]

            // compute upper intersection
            // Get costs/distances
            typedef boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
            edge_descriptor p = boost::edge(p1.node_id, p2.node_id, *g).first;
            edge_descriptor q = boost::edge(q1.node_id, q2.node_id, *g).first;

            float dist_p2_q1 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j) : RoadNetworkPoint::dist(p2,q1);
            float dist_p2_q2 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j + 1) : RoadNetworkPoint::dist(p2,q2);
            float cost_p = boost::get(boost::edge_weight, *g, p);
            float cost_q = boost::get(boost::edge_weight, *g, q);

            // Define Lines
            float l2_slope = cost_q / cost_p;
            float l2_intercept = (dist_p2_q1 + cost_p * (cell.i + 1) - cost_q * cell.j - epsilon)/cost_p;

            float l3_slope = -cost_q / cost_p;
            float l3_intercept = (dist_p2_q2 + cost_p * (cell.i + 1) + cost_q * (cell.j + 1) - epsilon)/cost_p;

            float intersect_x = (l3_intercept - l2_intercept) / (l2_slope - l3_slope);
            float intersect_y = l2_slope * intersect_x + l2_intercept;

            auto upper_right = (intersect_y >= cell.i + 1) ? FreeSpaceCell::null_pair : std::pair<float,float>(intersect_y, cell.i + 1);
            cell.setRight(cell.getRight().first, upper_right);

        }

    }

}

/*
    Computes given Free Space cell side
 */
void RoadNetworkPoint::computeFreeSpaceSide(RoadNetworkPoint& p1, RoadNetworkPoint& q1, RoadNetworkPoint& p2, RoadNetworkPoint& q2, FreeSpaceCell& cell, float epsilon, side s, std::shared_ptr<DistCache<RoadNetworkPoint>> cache){

    /*
        See createFreeSpace function for more info on the computational aspects
     */

    // Edge costs
    typedef boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
    edge_descriptor p = boost::edge(p1.node_id, p2.node_id, *g).first;
    edge_descriptor q = boost::edge(q1.node_id, q2.node_id, *g).first;

    float cost_p = boost::get(boost::edge_weight, *g, p);
    float cost_q = boost::get(boost::edge_weight, *g, q);

    if(s == BOTTOM){
        float dist_p1_q1 = (cache) ? cache->cacheDistance(cell.i, cell.j) : RoadNetworkPoint::dist(p1,q1);
        float dist_p1_q2 = (cache) ? cache->cacheDistance(cell.i, cell.j + 1) : RoadNetworkPoint::dist(p1,q2);
        float l1_slope = -cost_q / cost_p;
        float l1_intercept = -(dist_p1_q1 - cost_p * cell.i - cost_q * cell.j - epsilon)/cost_p;
        float l4_slope = cost_q / cost_p;
        float l4_intercept = -(dist_p1_q2 - cost_p * cell.i + cost_q * (cell.j + 1) - epsilon)/cost_p;

        auto b = line_line_intersection(l1_slope, l1_intercept, l4_slope, l4_intercept, cell.i, cell.j, false, false);
        cell.setBottom(b.first, b.second);
        return;
    }

    else if(s == LEFT){
        float dist_p1_q1 = (cache) ? cache->cacheDistance(cell.i, cell.j) : RoadNetworkPoint::dist(p1,q1);
        float dist_p2_q1 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j) : RoadNetworkPoint::dist(p2,q1);
        float l1_slope = -cost_q / cost_p;
        float l1_intercept = -(dist_p1_q1 - cost_p * cell.i - cost_q * cell.j - epsilon)/cost_p;
        float l2_slope = cost_q / cost_p;
        float l2_intercept = (dist_p2_q1 + cost_p * (cell.i + 1) - cost_q * cell.j - epsilon)/cost_p;

        auto l = line_line_intersection(l1_slope, l1_intercept, l2_slope, l2_intercept, cell.i, cell.j, true, false);
        return;
    }

    else if(s == RIGHT){
        float dist_p1_q2 = (cache) ? cache->cacheDistance(cell.i, cell.j + 1) : RoadNetworkPoint::dist(p1,q2);
        float dist_p2_q2 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j + 1) : RoadNetworkPoint::dist(p2,q2);
        float l3_slope = -cost_q / cost_p;
        float l3_intercept = (dist_p2_q2 + cost_p * (cell.i + 1) + cost_q * (cell.j + 1) - epsilon)/cost_p;
        float l4_slope = cost_q / cost_p;
        float l4_intercept = -(dist_p1_q2 - cost_p * cell.i + cost_q * (cell.j + 1) - epsilon)/cost_p;

        auto r = line_line_intersection(l3_slope, l3_intercept, l4_slope, l4_intercept, cell.i, cell.j, true, true);
        return;
    }

    else if(s == TOP){
        float dist_p2_q1 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j) : RoadNetworkPoint::dist(p2,q1);
        float dist_p2_q2 = (cache) ? cache->cacheDistance(cell.i + 1, cell.j + 1) : RoadNetworkPoint::dist(p2,q2);
        float l2_slope = cost_q / cost_p;
        float l2_intercept = (dist_p2_q1 + cost_p * (cell.i + 1) - cost_q * cell.j - epsilon)/cost_p;
        float l3_slope = -cost_q / cost_p;
        float l3_intercept = (dist_p2_q2 + cost_p * (cell.i + 1) + cost_q * (cell.j + 1) - epsilon)/cost_p;

        auto t = line_line_intersection(l2_slope, l2_intercept, l3_slope, l3_intercept, cell.i, cell.j, false, true);
        return;
    }

}