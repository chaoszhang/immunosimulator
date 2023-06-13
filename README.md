# immunosimulator
Dynamic IMmuno-SIMulator

## simulator.hpp
Generic simulator under the birth-death-mutation process.

### class Taxon
Generic class for an unique group of individuals.
Specimens in the same taxon share the same properties (eg. birth rate, occupancy, etc).
This class should not be modified. Instead, its derivative classes should be created.

### class ExampleTaxon
A minimum example of how to define a taxon.

### class Simulator
Generic simulator under the birth-death-mutation process with user defined taxon class.
This class and its dependents should not be modified and can be reused as is.

## main.cpp
Our adoption of the simulator for simulating the somatic hypermutation process.
The binary file produced by main.cpp depends on kmerFreq.txt and a parameter file produced by setup.py.

### class BCell
Our defination of Taxon in the context of BCR coding sequence evolution.

## kmerFreq.txt
K-mer Frequences for the 6-mer model.
Columns represent the frequences that the central nucleotide change to A, C, G, and T, respectively.

## setup.py
The script to produce the parameters for various simulation conditions.

## stat
Various scripts for benchmarking the performance of reconstruction tools using various metrics.
