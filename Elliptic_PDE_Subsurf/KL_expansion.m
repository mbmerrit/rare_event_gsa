function [KL] = KL_expansion(parameters, stochastic_dim)

log_kappa_ave    = parameters.log_kappa_ave; 
log_kappa_stddev = parameters.log_kappa_stddev;
nx = parameters.nx;
ny = parameters.ny;

w = WeightsCal(nx, ny);

%sum(w(:))
%--- stochastic distribution of kapa ---%
% build the correlation struct
%
corr.name = 'exp';          % the type of correlation function
cx = parameters.corrx;      % correlation length in x direction
cy = parameters.corry;      % correlation length in y direction 
corr.c0 = [cx cy];          % correlation length vector 
corr.sigma = log_kappa_stddev.^2;       % set a constant variance 

% compute the eigenvalues and eigenfunctions
quadrature.nodes = parameters.mesh;
quadrature.weights = w(:); 
trunc = stochastic_dim;
[KL] = randomfield_modified(corr, log_kappa_ave, quadrature, trunc); 

% compute the average standard deviation (should equal log_kappa_stddev) 
avgstd = sqrt(mean(KL.basis.^2 * KL.sv.^2));
%fprintf('Average std deviation of the field over the domain = %8.4f (exact stddev = %8.4f)\n', avgstd, log_kappa_stddev);


