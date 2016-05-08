# SF2568 Parallel Computations for Large Scale Systems
## Project - Spring 2016

Finding the optimal solution to the k-means clustering problem in general Euclidean space is NP-hard, 
and also NP-hard for a general number of clusters k (even in the plane). However, running the k-means algorithm can be
relatively fast so it is not uncommon to run it multiple times with different initial number of clusters and/or cluster means. In this
report we do not aim to find the optimal solution, but rather to explore how great a speedup can be achieved by making a parallel
implementation of the k-means algorithm.

This repository has both the serial version of the algorithm along with the parallel version written in C. In the folder pythonScripts is a script to generate the datasets.

To build the source code you can use the following command. This will build the source code and add it to a bin folder.
```
make
```
To clean the build code
```
make clean
```

