% Computes Sobol' indices for the 2D subsurface flow problem using Saltelli
% sampling algorithm (MC) and Subset simulation (SS) for rare event
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
%
% some initialization
%
N_SS = 5e2; % number of inner QoI evals for SS for each outer QoI eval
N_Sobol = 1e3; % number of Sobol samples in each dimension
ndim = 3; % dimension of SLSA param vector
tau_RE = 4.5; % rare event threshold
p0 = 0.1; % quantile for SS

% nominal parameter values
k1 = 0.4; % corr length x
k2 = 0.4; % corr length y
k3 = 0.8; % sigma
mu = [k1; k2; k3];   % mean of the parameters

%
% 10% pertrubation around mean
%
perturb = 0.1;
a = mu - perturb*mu;
b = mu + perturb*mu;

%
% parameter samples (generate U[-1,1] realizations)
%
XA = -1 + 2*rand(N_Sobol, ndim); 
XB = -1 + 2*rand(N_Sobol, ndim); 
A = repmat(0.5*(a+b)',N_Sobol,1) + 0.5*(b-a)'.*XA; % hyperparameters
B = repmat(0.5*(a+b)',N_Sobol,1) + 0.5*(b-a)'.*XB; % hyperparameters
qC = zeros(N_Sobol,ndim);

%
% sampling loop
%
for i = 1:N_Sobol
    [KLa, pde_data] = build_kle_fromparams(50, 50, A(i,1), A(i,2));
    [qA(i), q_eval, tau, theta] = SubSim(tau_RE, N_SS, p0, KLa, pde_data, A(i,3));
    
    [KLb, pde_data] = build_kle_fromparams(50, 50, B(i,1), B(i,2));
    [qB(i), q_eval, tau, theta] = SubSim(tau_RE, N_SS, p0, KLb, pde_data, B(i,3));
    disp(['finished QoI evaluation ', num2str(i), ' of ', num2str(N_Sobol)])
end
disp('done with A and B')

for k = 1:ndim % this substitutes one column of B into each A, then evaluates
   C = A;
   C(:,k) = B(:,k);  % here you are varying one parameter at a time
   for i = 1:N_Sobol
       [KLc, pde_data] = build_kle_fromparams(50, 50, C(i,1), C(i,2));
       [qC(i,k), q_eval, tau, theta] = SubSim(tau_RE, N_SS, p0, KLc, pde_data, C(i,3));
       disp(['finished QoI evaluation ', num2str(i), ' of ', num2str(N_Sobol)])
   end
   disp(['done with C',num2str(k)])
end

% Here is where we compute the Sobol' Indices - just first order and total
[SobolIndices] = get_sobol_indices(qA', qB', qC);

