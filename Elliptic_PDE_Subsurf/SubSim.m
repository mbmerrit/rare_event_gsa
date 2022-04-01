function [prob_re, q_eval, c, theta, counter] = SubSim(c_RE, N, p0, varargin)
% Subset Simulation code that evaluates the "compute_qoi" function. Also
% uses the modified Metropolis implementation given by MCMC_samp. 
% INPUT: the rare event threshold c_RE, number of samples per event region,
% and the fixed cond. prob. value for each region. It has been shown that
% the optimal p0 lies in [0.1, 0.3], see Sehic 2020.
% OUTPUT: rare event prob, function evaluations, failure levels for each 
% rare event region, and input Markov Chains.
%
% Implemented April 2021, Mike Merritt.
%

% Initialization
ndim = 126; % this is fixed by the user
Nc = ceil(N/(N*p0)); % length of an individual Markov Chain
q_eval = cell(1); % this cell array will grow for each region
q_eval_mat = cell(1); % same as above
theta_seeds = cell(1); % these are QoI inputs, as seeds for a MC at each level
counter = zeros(1, ceil(p0*N)); % counts the accepeted/rejected MCMC steps, or all QoI evals
gamma = 0.9; % this will be fixed for now, can also be chosen adaptively!
i = 1; % the failure region counter
if nargin > 5
    KLa = varargin{1}; % If you provide a full KL struct, we can use it
    pde_data = varargin{2}; % also need the PDE mesh data
    sigma = varargin{3};
end

% Generate initial samples, evaluate and sort them
counter(1,1) = N;
theta{i,i} = randn(N, ndim); % this distribution is problem specific
q_eval{i,i} = compute_qoi(theta{i,i}, 'fix_kle', KLa, pde_data, sigma)';
q_eval_mat{i} = q_eval{i,i};
% Define c1 as the p0 quantile and first event region
c(i) = quantile(q_eval_mat{i},1-p0); % the first threshold
theta_seeds{i} = theta{i,i}(q_eval_mat{i} > c(i), :); % seeds of next level of Markov Chains

% Main loop, uses MCMC, estimates new c and collects Markov Chain seeds
while c(i) < c_RE && i < 1e2  % chosen termination criterion 
    %disp(['level ', num2str(i), ' sampling done!'])
    %disp(['c threshold for level = ', num2str(c(i))])
    i = i + 1; % increment the levels only if you need to
% This approach is more complicated, but twice as fast from parallelization
    parfor j = 1:ceil(p0*N) % create N*p0 Markov Chains and store output
        [theta{i,j}, counter(i,j), q_eval{i,j}] ...
             = MCMC_samp(theta_seeds{i-1}(j,:), gamma, c(i-1), Nc, KLa, pde_data, sigma);
         theta_tmp{j} = theta{i,j}; % temporary storage for this level
         q_eval_tmp{j} = q_eval{i,j}; % temporary storage for this level 
    end
    theta_tmp_mat = vertcat(theta_tmp{:}); % temp storage, unzip cell array
    q_eval_mat{i} = vertcat(q_eval_tmp{:}); % same, just unzipping cells
    c(i) = quantile(q_eval_mat{i}, 1-p0);
    % then we collect the seeds of next level of Markov Chains
    theta_seeds{i} = theta_tmp_mat(q_eval_mat{i} >= c(i), :); 

% % This approach is simpler, but uses less parallelization and so is slower   
%     for j = 1:ceil(p0*N) % create N*p0 Markov Chains and store output
%          [theta{i}((j-1)*Nc+1:j*Nc,:), ~, q_eval{i}((j-1)*Nc+1:j*Nc)] ...
%              = MCMC_samp(theta_seeds{i-1}(j,:), gamma, c(i-1), Nc, KLa, pde_data);
%     end
%     c(i) = quantile(q_eval{i},1-p0);
%     theta_seeds{i} = theta{i}(q_eval{i} > c(i));
end
% Generate final samples
NF = sum(q_eval_mat{end} > c_RE);
% Post process, compute rare event prob
prob_re = p0^(i-1) * (NF/N);
q_eval = q_eval_mat; % output the matrix version for cleaner analysis