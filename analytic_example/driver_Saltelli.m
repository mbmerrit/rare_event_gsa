%
% This script implements the Saltelli MC approach for GSA with analytic
% expressions for computing the rare event probability. Cost of this method
% is dimension dependent. N controls the MC samples evaluated in each 
% dimension. The total cost is N(d + 2) for all main effect/total indices. 
% 
% Implemented April 2021, Mike Merritt. 
%
tic
%
% some initialization
%
if exist('N')
    ;
else
    N = 1e3; % number of samples for Sobol computation - in reality N(ndim + 2)
end
ndim = 10; % mean and variance for each RV
beta = 3; 

% nominal parameters
mu = [1; 2; 3; 4; 5];   % mean of the RVs
sigma = [10; 8; 6; 4; 2]; % variance of RVs
MU = [mu; sigma]; % mean of hyperparameters

%
% 10% pertrubation around mean
%
a = MU - 0.1*MU;
b = MU + 0.1*MU;

%
% parameter samples (generate U[-1,1] realizations)
%
%XA = -1 + 2*rand(N, ndim); 
%XB = -1 + 2*rand(N, ndim); 
XA = -1 + 2 * lhsdesign(N, ndim);  % sampling via LHS, can reduce variance
XB = -1 + 2 * lhsdesign(N, ndim);  % but drawing samples is more expensive

A = repmat(0.5*(a+b)',N,1) + 0.5*(b-a)'.*XA;
B = repmat(0.5*(a+b)',N,1) + 0.5*(b-a)'.*XB;
qC = zeros(N,ndim);

%
% sampling loop
%
[qA] = compute_qoi_analytic(A,beta);
fprintf('Finished sampling A / ');
[qB] = compute_qoi_analytic(B,beta);
fprintf('B / ');
for k = 1:ndim % this substitutes one column of B into each A, then evaluates
   C = A;
   C(:,k) = B(:,k);  % here you are varying one parameter at a time
   [ qC(:,k)] = compute_qoi_analytic(C,beta); 
   fprintf('C%d / ',k);
end

% Here is where we compute the Sobol' Indices - just first order and total;
% Or just accept an entire struct of Sobol index information
[SobolIndices] = get_sobol_indices(qA, qB, qC);
%SobolIndices

toc