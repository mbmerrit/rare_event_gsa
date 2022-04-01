function [qoi] = compute_qoi_analytic(X,tau)
N = size(X,1);

qoi = zeros(N,1);

parfor i = 1 : N


  % show progress every 100 samples
  if mod(i, 5e5)==0
     disp(i);
  end
  
  % This part of the QoI evaluation is hard coded for simplicity
  d = 5; % dimensions
  mu = X(i,1:d); sigma = X(i,d+1:2*d);
  mu_Z = -sum(mu) / sqrt(d); var_Z = sum(sigma) / d;
  P_RE = 1 - (1/2) * (1 + erf( ( tau - mu_Z) / sqrt(2 * var_Z)  ));
  qoi(i) = P_RE;

end