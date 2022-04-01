function [pce L] = uq_l1regression(pcdata, X, Y, tau, opts)
% purpose: compute PCE of a random variable using sparse regression
% input:   
%      pcdata --- struct containing basic PCE information 
%      X      --- sampleset of the input parameters  
%      Y      --- model evaluations at the sampling points 
%
% output: PC coefficients
 
%addpath('spgl1-1.9')

% get the needed structures
ndim = pcdata.ndim;
nord = pcdata.nord;
L = uq_get_forward_map(pcdata, X);
sigma = [];    % lasso option
x0  = [];
pce = spgl1(L, Y(:), tau, sigma, x0, opts);


