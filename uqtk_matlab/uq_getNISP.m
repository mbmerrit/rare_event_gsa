function K = uq_getNISP(pcdata, quadrature)
% function: uq_getNISP
%  
% purpose: Compute the NISP matrix which is used to 
%          project data into a PC basis.
%
% Usage: K = uq_getNISP(pcdata, quadrature)
%
% input:
%    pcdata     --- struct containing basic pc parameters    
%    quadrature --- stuct containing nodes and weights of multi-dimensional
%                   quadrature
%
%   
% output:
%    K      --- the NISP matrix
%       K_{ij} = W_j * Psi_i(X_j) / <Psi_i^2>
%         
% note: Given a vector F of evaluations of a function f
%       on quadrature nodes X_j, F_j = f(X_j), 
%       projection of X into a PC basis is accomplished 
%       through the matrix vector product
%            c = K * F
%       here c denotes the vector of spectral coeffs of f.
%
% Licensing:
%
%     Copyright (2011) Sandia Corporation. Under the terms of Contract
%     DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
%     retains certain rights in this software.
%
%     This file is part of The UQ Toolkit (UQTk)
%
%     UQTk is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     UQTk is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public License
%     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
%
%     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
%     Sandia National Laboratories, Livermore, CA, USA


% get the needed info from pcdata:
ndim = pcdata.ndim;
nord = pcdata.nord;
pc_type = pcdata.pc_type;
nPCTerms  = pcdata.nPCTerms;
psiMultiDSq = pcdata.psiMultiDSq;

% extract the mult-dimensional nodes and weights of quadrature
X = quadrature.nodes;
W = quadrature.weights;
nquad = quadrature.nquad;

% compute the NISP matrix
K = zeros(nPCTerms, nquad);
for j = 1 : nquad

    % multi-dimensional polynomials at quadrature nodes
    Psi_Xj = uq_PCBasis(pcdata, X(j, :));

    % jth row of K
    K(:, j) = W(j) * Psi_Xj(:) ./ psiMultiDSq(:);
end

