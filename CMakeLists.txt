cmake_minimum_required(VERSION 3.27)
project(traj_sim)

set(CMAKE_CXX_STANDARD 17)
find_package(GEOS CONFIG REQUIRED)

set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED)

if(Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})

endif()

add_custom_target(phl ALL COMMAND make run
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/External/pruned-highway-labeling
)

add_executable(traj_sim main.cpp
        External/pruned-highway-labeling/src/pruned_highway_labeling.cpp
        TrajModels/FreeSpaceCell.cpp
        TrajModels/EuclideanPoint2D.cpp
        TrajModels/RoadNetworkPoint.cpp
        TrajModels/EuclideanPoint3D.cpp
        DataLoading/DataLoaderEuclidean.cpp
        Utils/Utilities.cpp
        Experiments/AbstractMeasure.hpp
        DataLoading/DataLoaderRoadNetwork.cpp
        Experiments/BaselineAbstract.hpp
        ModifiedSimMeasures/ModifiedFrechet.hpp
        ModifiedSimMeasures/ModifiedHausdorff.hpp
        ModifiedSimMeasures/HausdorffIndexManager.hpp
        ModifiedSimMeasures/DTWPathManager.hpp
        ModifiedSimMeasures/ModifiedDTW.hpp
        ModifiedSimMeasures/CombinedSimilarityMeasures.hpp
        Experiments/CombinedMeasuresExperiments.hpp
        Utils/Bitmask.cpp
        Utils/Bitmask.hpp
        ModifiedSimMeasures/ExperimentalVersions.hpp
        ModifiedSimMeasures/ExperimentalEpsArrayVersions.hpp
        Experiments/SimilarityAnalysis.cpp
        Experiments/SimilarityAnalysis.hpp
)

target_link_libraries(traj_sim PRIVATE GEOS::geos_c)
