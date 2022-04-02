% Computes Sobol' indices for the 2D subsurface flow problem using Gauss 
% quadrature for building the PCE and Subset simulation (SS) for rare event
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
% Rare event initialization for subset simulation code
N_SS = 250; % number of inner QoI evals for SS for each outer QoI eval
tau_RE = 4.5;
p0 = 0.1;

% Nominal parameters for SLSA
k1 = 0.4; % corr length x
k2 = 0.4; % corr length y
k3 = 0.8; % sigma
mu = [k1; k2; k3];   % mean of the parameters
perturb = 0.1; % perturbation of parameters about the mean

% PCE quadrature setup
numberOf1DNodes = 7;
[quadrature] = uq_quadrature(ndim, numberOf1DNodes, pc_type);
Nq = quadrature.nquad;
X = quadrature.nodes;
W = quadrature.weights;
K = uq_getNISP(pcdata, quadrature); % transforms QoI evals into PCE coefs

X_transf = X .* mu' * perturb + mu'; % the input transformed to non standard intervals
int_length = 2*perturb*mu; % size of transformed interval for Uniform RVs

%
% sampling loop
%
QoI_evals = zeros(Nq, 1);
for i = 1 : Nq
    [KLa, pde_data] = build_kle_fromparams(50, 50, X_transf(i,1), X_transf(i,2));
    [QoI_evals(i), q_eval, tau, theta] = SubSim(tau_RE, N_SS, p0, KLa, pde_data, X_transf(i,3));
    disp(['finished QoI evaluation ', num2str(i), ' of ', num2str(Nq)])
end

pce_quad = K * QoI_evals;
toc

% This code can be used to sample the constructed PCE for verification
% N_ver = 1e5;
% X_ver = -1 + 2 * rand(N_ver, ndim);
% X_ver_tr = X_ver .* mu' * perturb + mu';
% 
% for i = 1 : N_ver
%    Y_pce_quad(i) = uq_evalpce(pcdata, pce_quad, X_ver(i,:));
% end

% Sobol index computation step
for i = 1 : ndim

   [T_quad(i)] = uq_getTotalSensitivity(pcdata, pce_quad, i);

end



