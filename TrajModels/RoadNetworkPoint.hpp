#pragma once
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "../External/pruned-highway-labeling/src/pruned_highway_labeling.h"

#include "../DistCache/DistCache.hpp"
#include "FreeSpaceCell.hpp"

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::property<boost::edge_weight_t, float>> Graph;
typedef std::pair< int, int > Edge;

/*
    Models points on a road network graph
 */
class RoadNetworkPoint {
private:

    static std::pair<fpair,fpair> line_line_intersection(float slope1, float intercept1, float slope2, float intercept2, int i, int j, bool vertical, bool upper);

public:

    static const int fs_intervals_per_side = 2;


    static std::shared_ptr<Graph> g;
    static std::shared_ptr<PrunedHighwayLabeling> phl;
    int node_id;

    explicit RoadNetworkPoint(int node_id ) : node_id(node_id) {};

    std::string to_string();

    static void setGraph(std::shared_ptr<Graph> graph){
        g = std::move(graph);
    };

    static void setPHL(std::shared_ptr<PrunedHighwayLabeling> pruned_highway_labeling){
        phl = std::move(pruned_highway_labeling);
    };

    static float dist(RoadNetworkPoint& i, RoadNetworkPoint& j);

    static FreeSpaceCell createFreeSpace(RoadNetworkPoint& p1, RoadNetworkPoint& q1, RoadNetworkPoint& p2, RoadNetworkPoint& q2, float epsilon, int i, int j, std::shared_ptr<DistCache<RoadNetworkPoint>> cache = NULL);
    static void propagateFreeSpace(RoadNetworkPoint& p1, RoadNetworkPoint& q1, RoadNetworkPoint& p2, RoadNetworkPoint& q2, FreeSpaceCell& cell, float epsilon, std::shared_ptr<DistCache<RoadNetworkPoint>> cache = NULL);
    static void computeFreeSpaceSide(RoadNetworkPoint& p1, RoadNetworkPoint& q1, RoadNetworkPoint& p2, RoadNetworkPoint& q2, FreeSpaceCell& cell, float epsilon, side s, std::shared_ptr<DistCache<RoadNetworkPoint>> cache = NULL);

    // Only declared, should remain unnused!
    static float dist_to_edge(RoadNetworkPoint& p1, RoadNetworkPoint& p2, RoadNetworkPoint& p3){return -1;}
    static float dist_of_bisector(RoadNetworkPoint& p1, RoadNetworkPoint& p2, RoadNetworkPoint& p3, RoadNetworkPoint& p4){return -1;}

};

