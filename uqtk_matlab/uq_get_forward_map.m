function K = uq_get_forward_map(pcdata, sampleset)
% function: uq_get_forward_map
%  
% Usage: K = uq_get_forward_map(pcdata, quadrature)
%
% input:
%    pcdata     --- struct containing basic pc parameters    
%    quadrature --- stuct containing nodes and weights of multi-dimensional
%                   quadrature
%
%   
% output:
%    K      --- the coefficient vecotor to data operator:
%       K_{jk} = Psi_k(X_j) 
%    
%    Note: for a pce vector c,  (K c)_j = sum_k  K_{jk} c_k = sum_k c_k Psi_k(X_j) = PCE(X_j).

% get the needed info from pcdata:
ndim = pcdata.ndim;
nord = pcdata.nord;
pc_type = pcdata.pc_type;
nPCTerms  = pcdata.nPCTerms;
psiMultiDSq = pcdata.psiMultiDSq;

X = sampleset; 
nquad = length(X); 

% compute the operator 
K = zeros(nquad, nPCTerms);
for j = 1 : nquad

    % multi-dimensional polynomials at quadrature nodes
    Psi_Xj = uq_PCBasis(pcdata, X(j, :));

    % jth row of K
    K(j, :) = Psi_Xj(:)'; 
end

