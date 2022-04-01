# Global sensitivity analysis of rare event probabilities

## Repository overview
This repository contains MATLAB code implementing and demonstrating the numerical experiments conducted in the paper:

Merritt, Michael, Alen Alexanderian, and Pierre Gremaud. ["Global sensitivity analysis of rare event probabilities."](https://arxiv.org/abs/2110.13974) arXiv preprint arXiv:2110.13974 (2021).

The code from the paper is organized into three directories:
- The analytic example of a linear limit state function, where global sensitivity analysis (GSA) is performed on the hyperparameters defining the normal random variables. 
- An elliptic PDE modeling Darcy flow in a 2D medium with unknown permeability. This problem is found in the context of nuclear waste repositories, where a KL expansion is commonly used to simulate the random field. GSA is performed on the KL hyperparameters. 
- A collection is third party codes, implementing various tasks, including sparse regression, Gauss quadrature, and construction of polynomial chaos expansions (PCE). 

## Global sensitivity analysis
Variance-based global sensitivity analysis is an essential tool for characterizing the importance of uncertain parameters in a mathematical model. 
