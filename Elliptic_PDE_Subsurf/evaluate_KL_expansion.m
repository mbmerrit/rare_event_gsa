function y = evaluate_KL_expansion(KL, xi)
% Y = EVALUATE_KL_EXPANSION(KL, XI)
%
%   Evaluates a Karhunen-Loeve expansion, 
%
%        Y = mu + sum_j sqrt(lambda_j) xi_j U_j 
%   
%   The eigenvalues (in fact the singular values) and the eigenvectors, and the mean
%   are stored in the struct KL

mu = KL.mean;
U  = KL.basis;
sv = KL.sv;
y = mu(:) + U * (sv(:) .* xi);


