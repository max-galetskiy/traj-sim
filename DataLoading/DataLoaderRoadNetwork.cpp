#include "DataLoaderRoadNetwork.hpp"

#include <fstream>
#include <iostream>
#include "boost/tokenizer.hpp"

/*
        Reads two files:
            - "graph_filepath" refers to the road network graph file (which comes in the form of a list of weighted edges)
            - "trajectory_filepath" refers to a csv containing trajectories on this road network graph

 */
void DataLoaderRoadNetwork::readFile(std::string graph_filepath, std::string trajectory_filepath, std::string phl_infilepath, std::string phl_outfilepath, unsigned int traj_cardinality, bool graph_header_present, bool traj_header_present, unsigned int source_column, unsigned int target_column, unsigned int distance_column, unsigned int traj_column){

    Graph g = readGraph(graph_filepath, graph_header_present, source_column, target_column, distance_column);
    std::shared_ptr<PrunedHighwayLabeling> phl = readPHL(phl_infilepath, phl_outfilepath);
    RoadNetworkPoint::setGraph(std::make_shared<Graph>(g));
    RoadNetworkPoint::setPHL(phl);

    readTrajectories(trajectory_filepath, traj_cardinality, traj_header_present, traj_column);
}

/*
        Reads a road network graph (which comes in the form of a list of weighted edges)
        Parameters decide which columns contain information on the locations of the source node, target node and the distance of the edge
 */
Graph DataLoaderRoadNetwork::readGraph(std::string filepath, bool header_present, unsigned int source_column, unsigned int target_column, unsigned int distance_column) {

    std::ifstream file(filepath);

    if(file){

        // Initialize Graph Object
        std::vector<Edge> edges = {};
        std::vector<float> weights = {};
        unsigned int max_node_id;

        // Read File
        std::string line;

        // Skip header
        if(header_present){
            std::getline(file,line);
        }

        while(std::getline(file, line)){

            // Parse columns
            std::vector<std::string> columns;
            boost::tokenizer<boost::escaped_list_separator<char>> line_tk(line, boost::escaped_list_separator<char>('\\', ';', '\"'));

            for (boost::tokenizer<boost::escaped_list_separator<char>>::iterator i(line_tk.begin()); i!=line_tk.end();++i)
            {
                columns.push_back(*i);
            }

            // Parse Edge
            unsigned int u = std::stoi(columns[source_column]);
            unsigned int v = std::stoi(columns[target_column]);
            float weight = std::stof(columns[distance_column]);

            edges.emplace_back(u,v);
            weights.push_back(weight);

            if(u > max_node_id || v > max_node_id){
                max_node_id = std::max(u,v);
            }

        }

        file.close();
        return {edges.begin(), edges.end(), weights.begin(), max_node_id};

    }else{
        std::throw_with_nested(std::runtime_error("Could not open file: " + filepath));
    }

}

/*
        Runs the Pruned Highway Labeling algorithm to generate labels which will later be used to query shortest paths in the road network (see README for more information on PHL)
 */
std::shared_ptr<PrunedHighwayLabeling> DataLoaderRoadNetwork::readPHL(std::string input_filepath, std::string output_filepath){
    auto phl = std::make_shared<PrunedHighwayLabeling>();
    phl->ConstructLabel(input_filepath.data());
    phl->StoreLabel(output_filepath.data());
    //phl.LoadLabel(output_filepath.data());
    return phl;
}

/*
        Reads "traj_cardinality" number of trajectories from a given csv file.
        Trajectories are expected to come in the following form:
        "node1;node2;node3;..."
 */
void DataLoaderRoadNetwork::readTrajectories(std::string filepath, unsigned int traj_cardinality, bool header_present, unsigned int traj_column) {
    std::ifstream file(filepath);

    if(file){

        // Initialize Trajectory Container
        typedef std::shared_ptr<Trajectory<RoadNetworkPoint>> traj_ptr;
        trajectories.reset();
        trajectories_size = 0;
        trajectories = std::shared_ptr<traj_ptr[]>(new traj_ptr[traj_cardinality]);

        // Read File
        std::string line;
        unsigned int n = traj_cardinality;

        // Skip header
        if(header_present){
            std::getline(file,line);
        }

        while((std::getline(file, line)) && (n > 0)){

            // Parse CSV columns
            std::vector<std::string> columns;
            boost::tokenizer<boost::escaped_list_separator<char>> line_tk(line, boost::escaped_list_separator<char>('\\', ';', '\"'));

            for (boost::tokenizer<boost::escaped_list_separator<char>>::iterator i(line_tk.begin()); i!=line_tk.end();++i)
            {
                columns.push_back(*i);
            }

            // Parse Trajectory
            std::string traj_string = columns[traj_column];

            // Trajectory is in the form "node,node,node,..."
            // Parse individual points
            std::vector<std::shared_ptr<RoadNetworkPoint>> points;
            boost::tokenizer<boost::escaped_list_separator<char>> pnt_tk(traj_string, boost::escaped_list_separator<char>('\\', ',', '\"'));

            for (boost::tokenizer<boost::escaped_list_separator<char>>::iterator i(pnt_tk.begin()); i!=pnt_tk.end();++i) {
                points.push_back(std::make_shared<RoadNetworkPoint>(std::stoi(*i)));
            }

            // Add trajectory to list
            trajectories[traj_cardinality - n] = std::make_shared<Trajectory<RoadNetworkPoint>>(points);
            trajectories_size++;

            n--;
        }

        file.close();
        return;

    }else{
        std::throw_with_nested(std::runtime_error("Could not open file: " + filepath));
    }
}
