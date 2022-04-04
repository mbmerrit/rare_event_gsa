% Computes Sobol' indices for the 2D subsurface flow problem using sparse
% regression approach to building the PCE with Latin hypercube sampling for
% the KLE hyperparameters and Subset simulation (SS) for rare event
% simulation. Each evaluation of the rare event probability requires one to
% build a KL expansion of the random field, solve the Darcy flow PCE, solve
% the time-evolution ODE, then extract the final time value. As a result,
% THIS CODE IS VERY TIME CONSUMING TO RUN. IT IS RECOMMENDED THAT THIS CODE
% BE RUN, IN PARALLEL, ON AN HPC. For a fast-running demonstration, see the
% analytic_example directory, those scripts run in seconds. 
%
% OUTPUT: Sobol' indices of rare event prob., function evaluations
%
% Implemented April 2021, Mike Merritt.
%
clear

tic
%
% some initialization
%
addpath ../uqtk_matlab
ndim = 3; % dimension of SLSA param vector
% PCE data
n_ord = 5; % PCE total polynomial order
pc_type = 'LEGENDRE';
pcdata = uq_pcset(n_ord, ndim, pc_type);
N_lhs = 1e3; % number of Latin hypercube samples, then use sparse regression
% Rare event initialization for subset simulation code
N_SS = 1e3; % number of inner QoI evals for SS for each outer QoI eval
tau_RE = 4.5;
p0 = 0.1;

% Nominal parameters for SLSA
k1 = 0.4; % corr length x
k2 = 0.4; % corr length y
k3 = 0.8; % sigma
mu = [k1; k2; k3];   % mean of the parameters
perturb = 0.1; % perturbation of parameters about the mean

X_lhs = -1 + 2 * lhsdesign(N_lhs, ndim); % LHS in interval [-1, 1]
X_transf = X_lhs .* mu' * perturb + mu'; % the input transformed to non standard intervals
int_length = 2*perturb*mu; % size of transformed interval for Uniform RVs

%
% sampling loop
%
QoI_evals = zeros(N_lhs, 1);
for i = 1 : N_lhs
    [KLa, pde_data] = build_kle_fromparams(50, 50, X_transf(i,1), X_transf(i,2));
    [QoI_evals(i), q_eval, tau, theta] = SubSim(tau_RE, N_SS, p0, KLa, pde_data, X_transf(i,3));
    disp(['finished QoI evaluation ', num2str(i), ' of ', num2str(N_lhs)])
end

toc

% Next, we post process the data with sparse regression codes
addpath ../uqtk_matlab/spgl1-1.9 % path to sparse regression codes
if exist('verb')
    ;
else
    verb = 1; % number of Latin hypercube samples for building PCE
end
opts = spgSetParms('verbosity',verb); % can be used to silence output
tau = 5e-2; % regularization coefficient for sparse regression
[pce_spreg, L] = uq_l1regression(pcdata, X_transf, QoI_evals, tau, opts);


% This code can be used to sample the constructed PCE for verification
% N_ver = 1e5; MU = mu;
% X_ver = -1 + 2 * rand(N_ver, ndim);
% X_ver_tr = X_ver .* MU' * perturb + MU';
% 
% for i = 1 : N_ver
%    Y_pce_spreg(i) = uq_evalpce(pcdata, PCE_spreg(1,:), X_ver(i,:));
% end

% Sobol index computation step
for i = 1 : ndim

    [T_spreg(i)] = uq_getTotalSensitivity(pcdata, pce_spreg, i);

end
