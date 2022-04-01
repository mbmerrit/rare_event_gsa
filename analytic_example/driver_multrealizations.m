clear all
%
% This code will take a significant amount of runtime, recommend to run in
% parallel and alter sampling conditions for your needs. Before using, take
% a look at the driver codes herein ('driver_Saltelli' and 'driver_spreg').
%
% Implemented April 2021, Mike Merritt. 
%
%%
 % Saltelli sampling approach with multiple realizations. Rare event
 % probability is computed analytically, then the standard MC method 
 % (pick and freeze) is then used to estimate Sobol' indices. 
Num_realizations = 100; N = 1e2; % number of realizations and Sobol samples
Sobol_S_Jan = zeros(Num_realizations,10);
Sobol_S_Sal = zeros(Num_realizations,10);
Sobol_T_Sal = zeros(Num_realizations,10);
Sobol_T_Sob = zeros(Num_realizations,10);
Sobol_T_Hom = zeros(Num_realizations,10);

for i = 1:Num_realizations
    disp(['---------------- Realization ',num2str(i), '-----------------'])
    driver_Saltelli
    Sobol_S_Jan(i, :) = SobolIndices.S_Jansen;
    Sobol_S_Sal(i, :) = SobolIndices.S_Saltelli;
    Sobol_T_Sal(i, :) = SobolIndices.T_Saltelli;
    Sobol_T_Sob(i, :) = SobolIndices.T_Sobol;
    Sobol_T_Hom(i, :) = SobolIndices.T_Homma;    
end

%% 
% PCE approach with Gauss-quadrature used to compute PCE coefficients.
% Total Sobol' indices are then computed as a post-process. This section 
% only computes a single realization of the PCE and Sobol' indices, to 
% be used as a reference solution, for comparison with the sparse PCE code.

mu = [1 2 3 4 5]; % vector of means
Sigma = [10 8 6 4 2]; % variances, no correlations here

addpath ../uqtk_matlab
addpath ../uqtk_matlab/spgl1-1.9
% PCE data
ndim = 10;
n_ord = 3; % PCE total polynomial order
pc_type = 'LEGENDRE';
pcdata = uq_pcset(n_ord, ndim, pc_type);
% Rare event initialization for subset simulation code
N_SS = 2.5e3; % number of inner QoI evals for SS for each outer QoI eval
tau = 3;

% Nominal parameters for SLSA
MU = [mu, Sigma]'; % mean of hyperparameters
perturb = 0.1; % perturbation of parameters about the mean
a = MU - perturb*MU;
b = MU + perturb*MU;

% PCE quadrature setup
numberOf1DNodes = 3;
[quadrature] = uq_quadrature(ndim, numberOf1DNodes, pc_type);
Nq = quadrature.nquad;
X = quadrature.nodes;
W = quadrature.weights;
K = uq_getNISP(pcdata, quadrature); % transforms QoI evals into PCE coefs

X_transf = X .* MU' * perturb + MU'; % the input transformed to non standard intervals
int_length = 2*perturb*MU; % size of transformed interval for Uniform RVs

QoI_evals = zeros(Nq, 1);
QoI_evals = compute_qoi_analytic(1, X_transf, tau);
pce_quad = K * QoI_evals;

for k = 1 : ndim
    [T_quad(k)] = uq_getTotalSensitivity(pcdata, pce_quad, k);
end

%% 
 % Sparse PCE sampling approach with multiple realizations. PCE generated
 % by LHS sampling of hyperparameters and sparse regression. Rare event
 % probability is computed analytically. Total Sobol' indices computed. 
Num_realizations = 100; N_lhs = 1e3;
Sobol_T_spreg = zeros(Num_realizations,10); 
verb = 0; % silences output of sparse regression code

for i = 1:Num_realizations
    disp(['---------------- Realization ',num2str(i), '-----------------'])
    driver_spreg
    for k = 1 : ndim
        [T_spreg(k)] = uq_getTotalSensitivity(pcdata, pce_spreg, k);
    end
    Sobol_T_spreg(i, :) = T_spreg;
end

%% Then plotting commands, need a sufficient amount of data for this section
bar(T_quad); hold on
errorbar(mean(Sobol_T_Sal), std(Sobol_T_Sal), 'ko', 'linewidth', 2)
errorbar(mean(Sobol_T_spreg), std(Sobol_T_spreg), 'ro', 'linewidth', 2)
title('Total index results','fontsize',16); 
legend('PCE with quadrature','Saltelli sampling, 1 std. dev.','PCE with regression, 1 std. dev.','fontsize',12)
xticklabels({'\mu_1','\mu_2','\mu_3','\mu_4','\mu_5','\sigma_1^2','\sigma_2^2','\sigma_3^2','\sigma_4^2','\sigma_5^2'})
axis tight; set(gca, 'fontsize', 20);

% High-fidelity solution, if one would like to compare with other results
%Sobol_T_true = [.01363, .052226, 0.1157123, .2029966, .312413369, .1619339, .104714, 0.059197, .026424895, .007037757];
%Sobol_S_true = [.0114749, .045744, .1031, .18246, .28437599, .146774, .0938566, .0528798, .023672, .005887];


