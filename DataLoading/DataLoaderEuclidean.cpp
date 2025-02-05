#include "DataLoaderEuclidean.hpp"

#include <fstream>
#include <iostream>
#include "boost/tokenizer.hpp"

/*
    Reads a csv file and saves traj_cardinality number of trajectories
    traj_column refers to the column where its points are stores

    Assumes a trajectory to be stored in the following fashion:
    "[[long,lat],[long,lat],...]"
 */
void DataLoaderEuclidean::readFile(const std::string& filepath, unsigned int traj_cardinality, bool header_present, unsigned int traj_column) {

        std::ifstream file(filepath);

        if(file){

            // Initialize Trajectory Container
            typedef std::shared_ptr<Trajectory<EuclideanPoint2D>> traj_ptr;
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
                boost::tokenizer<boost::escaped_list_separator<char>> line_tk(line, boost::escaped_list_separator<char>());

                for (boost::tokenizer<boost::escaped_list_separator<char>>::iterator i(line_tk.begin()); i!=line_tk.end();++i)
                {
                    columns.push_back(*i);
                }


                // Parse Trajectory
                std::string traj_string = columns[traj_column];

                // Trajectory is in the form [[long,lat], [long,lat], ...]]
                // First we replace all instances of "[" by double quotes and also remove the beginning end ending characters
                traj_string.erase(0,1);
                traj_string.erase(traj_string.length() - 1, 1);
                std::replace(traj_string.begin(), traj_string.end(), '[', '"');
                std::replace(traj_string.begin(), traj_string.end(), ']', '"');

                // Then we parse the trajectory as if itself were a csv file
                std::vector<std::shared_ptr<EuclideanPoint2D>> points;
                boost::tokenizer<boost::escaped_list_separator<char>> pnt_tk(traj_string, boost::escaped_list_separator<char>());

                for (boost::tokenizer<boost::escaped_list_separator<char>>::iterator i(pnt_tk.begin()); i!=pnt_tk.end();++i)
                {
                    const std::string& curr_pnt = *i;  // has the form "long,lat"  (including the double quotes!)

                    unsigned long long comma_index = curr_pnt.find(',');

                    float longitude = std::stof(curr_pnt.substr(0, comma_index));
                    float latitude = std::stof(curr_pnt.substr(comma_index + 1));

                    auto p = std::make_shared<EuclideanPoint2D>(std::pair<float,float>(longitude, latitude), true);
                    points.push_back(p);

                }


                // Add trajectory to list
                trajectories[traj_cardinality - n] = std::make_shared<Trajectory<EuclideanPoint2D>>(points);
                trajectories_size++;

                n--;
            }

            file.close();
            return;

        }else{

            std::throw_with_nested(std::runtime_error("Could not open file: " + filepath));

        }

}
