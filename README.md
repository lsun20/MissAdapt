# MissAdapt

This repository contains the code to replicate all numerical results in our paper "Adapting to Misspecification".  It also provides MATLAB code implementing the adaptive estimator, its soft-thresholding approximation, and their risk limited variants that we propose in the paper.

## Introduction
The outputs and the plots of the numerical examples are included in their respective folders, along with replication code for the original result.

All Matlab scripts for implementing the adaptive estimator are included in the folder `Matlab/`.   

## Matlab scripts

`/lookup_tables/` contains pre-tabulated solutions to the adaptation problem over a grid of correlation coefficients.

`adaptive_estimate.m` : a function that takes the input of the original results and outputs the adaptive estimates via an instantaneous lookup table.

`example.m` : an example script that calls `adaptive_estimate.m` and output the adaptive results in a formatted csv file.
 
