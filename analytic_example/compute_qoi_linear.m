function [qoi_tau] = compute_qoi_linear(X, hps)

    N = size(X,1);
    qoi_tau = zeros(1,N);
    % Initialize with user defined data
    
    % improves performance by switching between series and parallel
    if N == 1
        parforArg = 0;
    else
        parforArg = Inf;
    end
    
ndim = length(hps);    
mu = hps(1:ndim/2); sigma = hps(ndim/2+1:ndim); % hyperparameters
qoi_tau = -sum(X .* sqrt(sigma) + mu,2) / sqrt(ndim/2); % linear example QoI
  
end