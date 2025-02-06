# Trajectory Similarity Measures - Computation and Baselines in C++
When delving into the research field of moving object databases and specifically trajectory data itself, one is sooner or later bound to come into contact with trajectory similarity measures. The study of those is particularly promising, because both the qualitative and efficiency-related aspects of such similarity measures are far from being easily analyzable.

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
The following capabilities are provided in the project under these directories: 

### TrajModels
Here, an abstract trajectory model and classes for specific types of trajectory points (i.e. Euclidean and Road Network points) are provided. Additionally, functions computing the Free Space datastructure used in the classic Alt et Godau algorithm for the Fréchet distance are found here.

### SimMeasures
Here, functions for similarity measure computation and decision problem solving (i.e. checking whether Measure(T1,T2) <= Epsilon) are provided. For each measure, care was taken in implementing the most efficient computation algorithms found in current literature. In particular, the following measures are provided: Fréchet, discrete Fréchet, Hausdorff, DTW, LCSS, EDR, ERP.

### DataLoading
Basic loading of either Euclidean or Road Network trajectories from a .csv file (examples for input data can be seen in the `/Datasets` directory).

### DistCache
Container to cache point-wise distances of a pair of trajectories. This is useful if multiple similarity measures are to be computed on the same pair of trajectories (especially for the case of Road Network trajectories).

### Experiments
Abstract classes support the evaluation/benchmarking of similarity measures on a larger dataset. Additionally, the experimental setup from the evaluation of my Bachelor's thesis can be found here.

## Pointers to Literature
Since this I am publishing this repo to assist people who are working or wanting to work with trajectory similarity measures, I also felt inclined to include this essential list of publications that revolves around the classic similarity measures:
* [Fréchet Distance](https://doi.org/10.1142/S0218195995000064)
* [Discrete Fréchet Distance](https://www.researchgate.net/profile/Thomas-Eiter-2/publication/228723178_Computing_Discrete_Frechet_Distance/links/5714d93908aebda86c0d1a7b/Computing-Discrete-Frechet-Distance.pdf)
* [Hausdorff Distance](https://doi.org/10.1109/34.232073)
* [Dynamic Time Warping Distance](https://dl.acm.org/doi/10.5555/3000850.3000887)
* [Edit Distance on Real Sequences](https://doi.org/10.1145/1066157.1066213)
* [Edit Distance with Real Penalty](https://doi.org/10.5555/1316689.1316758)
* [Longest Common Subsequence](10.1109/ICDE.2002.994784)

Additionally, also see:
* [Adaption of the Fréchet Distance to Road Networks](https://doi.org/10.1007/978-3-642-24983-9_7)
* [For a more efficient computation of the Hausdorff Distance](https://doi.org/10.1109/TPAMI.2015.2408351)

## Authors and Acknowledgement
This project was written in its entirety by Max Galetskiy with the exception of the functions provided by Boost, GEOS or the PHL implementation as explained in the "Dependencies" section of this documentation.

The Bachelor's thesis from which this project originated from, was supervised by [Theodoros Chondrogiannis](https://dbis.uni-konstanz.de/people/people/researchers/chondrogiannis/) and [Johann Bornholdt](https://dbis.uni-konstanz.de/people/people/researchers/bornholdt/) of the [Database and Information Systems group](https://dbis.uni-konstanz.de/) under [Michael Grossniklaus](https://dbis.uni-konstanz.de/people/people/grossniklaus/).
