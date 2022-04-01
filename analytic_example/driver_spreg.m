%
% This script implements the sparse regression PCE approach with analytic
% expressions for computing the rare event probability. SPGL1 is used for
% sparse regression. No subset simulation here and 'N_lhs' controls the
% number of Latin hypercube samples used to build the PCE. 
% 
% Implemented April 2021, Mike Merritt. 
%
tic

mu = [1 2 3 4 5]; % vector of means
Sigma = [10 8 6 4 2]; % variances, no correlations here
addpath ../uqtk_matlab % path to UQ and quadrature codes
addpath ../uqtk_matlab/spgl1-1.9 % path to sparse regression codes
% PCE data
%load('PCE_data.mat')
ndim = 10;
n_ord = 2; % PCE total polynomial order
pc_type = 'LEGENDRE';
pcdata = uq_pcset(n_ord, ndim, pc_type);
% Rare event initialization
c_RE = 3;

% Nominal parameters for SLSA
MU = [mu, Sigma]'; % mean of hyperparameters
perturb = 0.1; % perturbation of parameters about the mean
a = MU - perturb*MU;
b = MU + perturb*MU;

% PCE LHS setup
if exist('N_lhs')
    ;
else
    N_lhs = 1e3; % number of Latin hypercube samples for building PCE
end
X_lhs = -1 + 2 * lhsdesign(N_lhs, ndim); % LHS in interval [-1, 1]
X_transf = X_lhs .* MU' * perturb + MU'; % the input transformed to non standard intervals

QoI_evals_lhs = zeros(N_lhs, 1);
QoI_evals_lhs = compute_qoi_analytic(X_transf, c_RE);

toc

if exist('verb')
    ;
else
    verb = 1; % number of Latin hypercube samples for building PCE
end
opts = spgSetParms('verbosity',verb); % can be used to silence output
tau = 5e-2; % regularization coefficient for sparse regression
[pce_spreg, L] = uq_l1regression(pcdata, X_lhs, QoI_evals_lhs, tau, opts);

