function [KL] = randomfield_modified(corr, mean, quadrature, trunc)
%% This is a simplified version of the code by P. Constantine and Q. Wang, 
%% I have removed the extra things we have no need for. 
%
%RANDOMFIELD Generates realizations of a Gaussian random field.
%
% KL = randomfield_modified(corr, mean, quadrature)
%   
% Outputs:
%
%   KL:         A struct containing the components of a Karhunen-Loeve (KL)
%               expansion of the random field. See below for a more detailed
%               description.
%
% Required inputs:
%   corr:       A struct containing the correlation information. See below
%               for details of the corr struct.
%
%   mean:       A user supplied mean of the random field. Must match the
%               values of the given mesh. (Default zeros)
%     
%   quadrature: a struct that contains the following fields:
%               quadrature.nodes     the computational grid (mesh)
%               quadrature.weights   the weights (i.e. quadrature weights) 
%
%
%   trunc:      Truncation level for the KL expansion.
%
%
%
% The output struct 'KL' with components of the Karhunen-Loeve representation
% of the random field has the following fields.
%   KL.mean:    The mean of the random field. 
%
%   KL.bases:   The eigenvectors covariance matrix.
%
%   KL.sv:      The square root of the eigenvalues of the covariance
%               matrix. Plot these on a semilog scale to examine their
%               decay.
%
% The input struct 'corr' may contain the following fields.
%   corr.name:  Specifies the correlation type from 'gauss', 'exp', or
%               'turbulent'.
%
%   corr.c0:    The scaling parameters for the correlation function. c0 may
%               be a scalar for isotropic correlation or a vector for
%               anisotropic correlation. In the anisotropic case, the
%               vector must have d elements, where d is the dimesion of a
%               mesh point.
%
%   corr.c1:    The second scaling parameters for the 'turbulent' 
%               correlation function. Not used with 'gauss' or 'exp'.
%
%   corr.sigma: The variance scaling parameter. May be a scalar or a vector
%               the size of the mesh. 
%
%
% References:
%   M. Davis, "Production of Conditional Simulations via the LU Triangular
%       Decomposition of the Covariance Matrix". Mathematical Geology,
%       1987.
%
% Copyright 2011 Qiqi Wang (qiqi@mit.edu) and Paul G. Constantine 
% (paul.constantine@stanford.edu).

if nargin<2, error('Not enough inputs.'); end

mesh = quadrature.nodes;
w = quadrature.weights;
nx=size(mesh,1);

% set default values.
nsamples=1;
data=[];
filter=1;
spthresh=0.0001;
X=[];
lowmem = 0;

% error checking the correlation struct
if ~isfield(corr,'name'), error('corr.name must be gauss, exp, or turbulent'); end
if ~isfield(corr,'c0'), error('corr.c0 missing.'); end
if ~isfield(corr,'c1')
    if strcmp(corr.name,'turbulent') 
        error('corr.c1 missing for turbulent correlation.'); 
    else
        corr.c1=[];
    end
end

% compute the correlation matrix
C = correlation_fun(corr, mesh, [], spthresh);

% compute the transform
%opts.issym=1;
opts.isreal=1;
opts.maxit=max(5*trunc,300);
opts.disp=0;
opts.v0=ones(nx,1);

W_sqrt = diag(sqrt(w));
K = W_sqrt * (C * W_sqrt);

K = 0.5 * (K + K'); % ensures K is symmetric, due to numerical loss of symmetry
[U, S, eigflag] = eigs(K,trunc,'lm',opts);

if eigflag==0
    fprintf(' eigs converged!\n');
else
    fprintf(' eigs convergence flag: %d\n',eigflag);
end

ev = abs(diag(S));
[ev, indz]=sort(ev, 'descend');
U = U(:, indz);

% weight eigenvectors:
U = W_sqrt \ U;

% construct the KL representation
KL.mean = mean;
KL.basis = U; 
KL.sv = sqrt(ev);
KL.C = C;
KL.nKL = trunc;
KL.K = K;
KL.trace = trace(K);
KL.r = cumsum(ev)./KL.trace;
