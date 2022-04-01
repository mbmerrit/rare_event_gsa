function [C, R, Ubar, Uc] = get_pressure_covariance(U)

% I assume rows of U are realizations of the process.

%
% compute the covariance matrix
%
Ns = size(U, 1);
disp('computing the covariance operator ...');
Ubar = mean(U);

%
% center the process
%
Uc = zeros(size(U));
for i = 1 : Ns
   Uc(i,:) = U(i, :) - Ubar;
end

% compute covariance matrix (this is the discretized covariance function)
C = cov(U);
R = corrcov(C);

