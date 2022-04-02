clear
%
% Exmaple using the simple linear limit state function with d independent
% normal random variables. The rare event probability includes the error
% function, but this is still a good, mostly analytic test. QoI: P_RE(mu, Sigma)
% 
n = 1e6; % number of realizations for each test, MAY NOT NEED
d = 5; % dimension, number of RVS in the sum
mu = [1 2 3 4 5]; % vector of means
sigma = [10 8 6 4 2]; % variances, no correlations here

%%
A = randn(n, d);
X = (A .* sqrt(sigma)) + mu; % transformed variables
Z = -sum(X, 2) / sqrt(d); % a new normally dist random variable
mu_Z = -sum(mu) / sqrt(d); var_Z = sum(sigma) / d;
disp(['Sum of means: ',num2str(mu_Z), ', and empirical sum: ', num2str(mean(Z))])
disp(['Variance of sum: ',num2str(var_Z), ', and empirical var: ', num2str(var(Z))])

% Given vectors of mus and sigmas, this computes the rare event prob
beta = 3; % rare event threshold
P_RE = 1 - (1/2) * (1 + erf( ( beta - mu_Z) / sqrt(2 * var_Z)  ));

[fi xi] = ksdensity(Z); plot(xi, fi, 'linewidth', 3); hold on
plot([beta, beta], [0 max(fi)], 'Linewidth', 2); title('PDF of QoI, q', 'fontsize',16)
disp(['Rare event probability = ', num2str(P_RE), ' and numerically = ', num2str(mean(Z > beta))])

%% First the rare event estimation step
close
N_SS = 100000;
HPs = [mu, sigma];
p0 = 0.1;
[P_RE, q_eval, c, theta] = SubSim_testprob(beta, N_SS, p0, HPs);
disp(['Subset simulation estimate of the rare event probability = ', num2str(P_RE)])

%% Now the Sobol step, compute sensitivity of P_RE to hyperparameters
% These should be considered the high-fidelity reference results. For a
% fair comparison of these codes with the PCE codes, reduce evaluations
N = 1e4; % with LHS sampling, this may be expensive, consider altering 
driver_Saltelli

%% Build a PCE for Prob_RE using Gauss quadrature, then do GSA from the PCE

tic
%
% some initialization
%
addpath ../uqtk_matlab
addpath ../uqtk_matlab/spgl1-1.9 % path to sparse regression codes
% PCE data
ndim = 10;
n_ord = 2; % PCE total polynomial order
pc_type = 'LEGENDRE';
pcdata = uq_pcset(n_ord, ndim, pc_type);
% Rare event initialization for subset simulation code
N_SS = 5e2; % number of inner QoI evals for SS for each outer QoI eval
p0 = 0.1;
c_RE = 3;

% Nominal parameters for SLSA
mu = [1 2 3 4 5]; % vector of means
sigma = [10 8 6 4 2]; % variances, no correlations here
MU = [mu, sigma]; % mean of hyperparameters
perturb = 0.1; % perturbation of parameters about the mean
a = MU - perturb*MU;
b = MU + perturb*MU;

% PCE quadrature setup
numberOf1DNodes = 2;
[quadrature] = uq_quadrature(ndim, numberOf1DNodes, pc_type);
Nq = quadrature.nquad;
X = quadrature.nodes;
W = quadrature.weights;
K = uq_getNISP(pcdata, quadrature); % transforms QoI evals into PCE coefs

X_transf = X .* MU * perturb + MU; % the input transformed to non standard intervals
int_length = 2*perturb*MU; % size of transformed interval for Uniform RVs

%
% sampling loop
%
QoI_evals_quad = zeros(Nq, 1);
for i = 1 : Nq
    if mod(i, 1e2)==0
        disp(['Function evaluation: ', num2str(i)])
    end
    % Need to reformulate the rare event problem BELOW the threshold
    [QoI_evals_quad(i), q_eval, c, theta] = SubSim_testprob(c_RE, N_SS, p0, X_transf(i,:));
end

pce_quad = K * QoI_evals_quad;
toc

for k = 1 : ndim
    [T_pce_quad(k)] = uq_getTotalSensitivity(pcdata, pce_quad, k);
end

% verification step - this blopck can sample the PCE for comparing PDFs

% N_ver = 1e5;
% X_ver = -1 + 2 * rand(N_ver, ndim);
% X_ver_tr = X_ver .* MU' * perturb + MU';
% 
% for i = 1 : N_ver
%    Y_pce_quad(i) = uq_evalpce(pcdata, pce_quad, X_ver(i,:));
% end

%% Latin hypercube sampling method with sparse regression for PCE coefs
%
tic
ndim = 10; % dimension of SLSA param vector
% PCE data
n_ord = 2; % PCE total polynomial order
pc_type = 'LEGENDRE';
pcdata = uq_pcset(n_ord, ndim, pc_type);
N_lhs = 1e3; % number of Latin hypercube samples, then use sparse regression
% Rare event initialization for subset simulation code
N_SS = 5e2; % number of inner QoI evals for SS for each outer QoI eval
c_RE = 3;
p0 = 0.1;

% Nominal parameters for SLSA
MU = [mu, sigma]; % mean of hyperparameters
perturb = 0.1; % perturbation of parameters about the mean
a = MU - perturb*MU;
b = MU + perturb*MU;

X_lhs = -1 + 2 * lhsdesign(N_lhs, ndim); % LHS in interval [-1, 1]
X_transf = X_lhs .* MU * perturb + MU; % the input transformed to non standard intervals
int_length = 2*perturb*MU; % size of transformed interval for Uniform RVs

%
% sampling loop
%
QoI_evals_lhs = zeros(N_lhs, 1);
for i = 1 : N_lhs
    if mod(i, 1e2)==0
        disp(['Function evaluation: ', num2str(i)])
    end
    % Need to reformulate the rare event problem BELOW the threshold
    [QoI_evals_lhs(i), q_eval, c, theta] = SubSim_testprob(c_RE, N_SS, p0, X_transf(i,:));
end
toc

verb = 1; % number of Latin hypercube samples for building PCE
opts = spgSetParms('verbosity',verb); % can be used to silence output
tau = 5e-2; % regularization coefficient for sparse regression
[pce_spreg, L] = uq_l1regression(pcdata, X_lhs, QoI_evals_lhs, tau, opts);

for k = 1 : ndim
    [T_spreg_lhs(k)] = uq_getTotalSensitivity(pcdata, pce_spreg, k);
end

%% Plot all results together
figure
bar([SobolIndices.T_Saltelli; T_pce_quad; T_spreg_lhs]')
title('Total index results','fontsize',16); 
legend('Saltelli sampling','PCE with quadrature','PCE with regression','fontsize',12)
xticklabels({'\mu_1','\mu_2','\mu_3','\mu_4','\mu_5','\sigma_1^2','\sigma_2^2','\sigma_3^2','\sigma_4^2','\sigma_5^2'})


