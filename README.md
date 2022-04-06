# Global sensitivity analysis of rare event probabilities

## Repository overview
This repository contains MATLAB code implementing and demonstrating the numerical experiments conducted in the paper:

Merritt, Michael, Alen Alexanderian, and Pierre Gremaud. ["Global sensitivity analysis of rare event probabilities."](https://arxiv.org/abs/2110.13974) arXiv preprint arXiv:2110.13974 (2021).

The code from the paper is organized into three directories:
- The analytic example of a linear limit state function, where Global Sensitivity Analysis (GSA) is performed on the hyperparameters defining the normal random variables. 
- An elliptic PDE modeling Darcy flow in a 2D medium with unknown permeability. This problem is found in the context of nuclear waste repositories, where a KL expansion is commonly used to simulate the random field. GSA is performed on the KL hyperparameters. 
- A collection is third party codes, implementing various tasks, including sparse regression, Gauss quadrature, and construction of polynomial chaos expansions (PCE). 

## Global sensitivity analysis
Variance-based global sensitivity analysis is an essential tool for characterizing the importance of uncertain parameters in a mathematical model. The Sobol' index is a popular GSA metric, where the sensitivity of an uncertain parameter is described according to its relative contribution to the model variance. Another common method for performing GSA is the use a polynomial surrogate model, which are also known as Polynomial Chaos Expansions (PCE). Given a PCE surrogate for a model, Sobol' indices can be computed at a negligible cost, making the PCE very effective in the context of GSA. 

## Rare event simulation
A rare event probability (or failure probability) quantifies the probability of some extreme event or failure in a system. Common examples are structural failure, tsunami damage prediction, or failure prediction in a nuclear waste repository. In each of these scenarios, uncertain model quantities are at play and assumptions are typically made about the distributional characteristics of said quantities. In this context, while GSA is an attractive tool, it represents a major computational challenge. By their naturem, estimating the probability of a rare event is inherently challenging. Through the use of advanced rare event simulation techniques and PCE, we aim to efficiently perform GSA on rare event probabilities. 

For additional details, contact Michael Merritt at mbmerrit@ncsu.edu
