# Simulation study code for "A potential outcomes approach to defining and estimating gestational age-specific exposure effects during pregnancy"

This repository contains code to implement the methods described in the paper "A potential outcomes approach to defining and estimating gestational age-specific exposure effects during pregnancy" by Schnitzer et al. (2021).

`datagen_func.R` contains the data generation function used in the simulation study. It has an option to generate counterfactual data under a fixed treatment initiation time (assigned to all subjects). 

`estimators_func.R` contains four functions, corresponding to the implementation of the target trials versions of the estimators IPW, G-computation, and TMLE described in the paper. The fourth function is the "standard" IPW method that ignores delivery time which was used in the paper.

`single_dataset_example.md` is a markdown file that demonstrates how a single dataset can be generated and the four functions applied. The file also details the specifications of the data generation parameters used in the paper's simulation study.
