function [seed_new, counter, q_eval] = MCMC_samp(seed, gamma, c, num_steps, hps)
% Implementation of a modified Metrpolis-Hastings algorithm for sampling
% from a target distribution phi(*|F_\ell-1). 
% INPUT: Given a random vector of seeds, with dimension ndim, gamma is 
% the correlation parameter in [0,1], c is the rare event threshold
% that the target must meet, and we give the number of steps in the chain. 
% OUTPUT: a Markov Chain of accepted steps in the parameter space.
% NOTE: more generally, the correlation (gamma) can be a vector, which can 
% be chosen adaptively, componentwise or as a scalar, to avoid rejection

% Implemented by Mike Merritt, April 2021.

[MC_size, ndim] = size(seed); 
counter = 0; % binary vector tabulating success of search
q_eval = zeros(num_steps, 1); % stores the QoI evaluations for the MC

% generate new candidate samples for the Markov Chain
q_eval(1) = compute_qoi_linear(seed, hps); % do not need this step, already have it
for i = 1:num_steps
    seed_cand = gamma*seed(end,:) + (sqrt(1-gamma^2) * randn(1, ndim));
    q_eval(i+1) = compute_qoi_linear(seed_cand, hps);
    if q_eval(i+1) > c
        seed = [seed; seed_cand]; %disp(['accepted with q = ',num2str(q_eval(i+1))])
    else
        seed = [seed; seed(end,:)]; %disp(['rejected with q = ',num2str(q_eval(i+1))])
        q_eval(i+1) = q_eval(i);
    end
end
seed_new = seed(2:end,:); % don't keep the first seed, from previous MC
q_eval = q_eval(2:end); % same with QoI evals, don't keep the first
counter = num_steps;


