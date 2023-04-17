# MissAdapt

This repository contains the code to replicate all numerical results in our paper "Adapting to Misspecification".  It also provides MATLAB and R code implementing the adaptive estimator, its soft-thresholding approximation, and their risk limited variants that we propose in the paper.

## Introduction
The outputs and the plots of the numerical examples are included in their respective folders, along with replication code for the original result.

All Matlab scripts for implementing the adaptive estimator are included in the folder `Matlab/`.   

All R scripts for implementing the adaptive estimator are included in the folder `R/`.  

## Matlab scripts

`/lookup_tables/` contains pre-tabulated solutions to the adaptation problem over a grid of correlation coefficients.

`adaptive_estimate.m` : a function that takes the input of the pre-tabulated solutions and outputs the adaptive estimates via an instantaneous lookup table.

`minimax_locus_plot.m` : a function that takes the input of the pre-tabulated solutions and outputs the locus of B-minimax estimates via an instantaneous lookup table.  The oracle risk function (risk of each B-minimax estimator) as well as the risk function of the adaptive estimator are also plotted.

`example.m` : an example script that calls `adaptive_estimate.m` and output the adaptive results in a formatted csv file; it also calls `minimax_locus_plot.m` and output the B-minimax locus plot along with risk functions.

## R scripts

`example.R` : an example script that reads the pre-tabulated solutions and outputs the adaptive estimates via an instantaneous lookup table.
 
