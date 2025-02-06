# Trajectory Similarity Measures - Computation and Baselines in C++
When delving into the research field of moving object databases and specifically trajectory data itself, one is sooner or later bound to come into contact with trajectory similarity measure. The study of those is particularly promising, because both the qualitative and efficiency-related aspects of such similarity measures are far from being easily analyzable.

## About
This project implements the computation of classic similarity measures (i.e. Fréchet, discrete Fréchet, Hausdorff, DTW, LCSS, EDR, ERP) in C++ for both Euclidean and Road Network trajectories. To my knowledge, a library combining all of these measures, has not yet been open sourced which is why I decided to publish my implementation. Originally, this project was part of my bachelor's thesis "A Framework and Optimizations for Trajectory Similarity Predicates" at the [University of Constance](https://www.uni-konstanz.de/en/). Hence, this repository also contains some experimental and benchmarking-related code.

## Dependencies
This code depends on the following C/C++ libraries:
* [GEOS](https://libgeos.org/)
* [Boost](https://www.boost.org/)

Regarding the installation of those, please refer to the according documentation. For a setup on Windows, the [vcpkg package manager](https://vcpkg.io/en/) is highly recommended.

Additionally, this project also uses [this implementation of the Pruned Highway Labeling algorithm](https://github.com/kawatea/pruned-highway-labeling), but the source code for those is already integrated within the project since some minor modifications had to be introduced to update the project to the C++17 standard. Hence, the PHL implementation does not need to be built seperately.

## How to build
To use, this project must be built from source. Both CMake and a C++ compiler supporting the C++17 standard are required. If these requirements are statisfied, then the `CMakeList.txt` file handles the building of the project. 

## Modules and Capabilities


## Pointers to Literature


## Authors and Acknowledgement
