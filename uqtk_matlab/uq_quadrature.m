function [quadrature] = uq_quadrature(ndim, numberOf1DNodes, pc_type)
% function: uq_quadrature.m
% purpose: returns nodes and weights of multi-
%          dimensional quadrature to be used when doing NISP.
%          Letting ndim be the stochastic dimension, this routine
%          returns the nodes and weights of a ndim-dimensional
%          fully tensorized Gauss quadrature.
%
%          This routine can be updated for future extensions 
%          enabling sparse / adaptive quadrature for NISP.
%
% usage: quadrature = uq_quadrature(ndim, numberOf1DNodes, pc_type)
%
% input:
%    ndim           :   stochastic dimension
%    numberOf1DNodes:   number of nodes in each direction
%    pc_type        :   the type of the PC basis 
%
% output:
%
%   quadrature      :   struct containing the quadrature data
%      quadrature.nodes     
%      quadrature.weights    
%      quadrature.nquad
% 
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



% lookup table for the quadrature points
ind = uq_quadtable(numberOf1DNodes, ndim);
[nquad unused] = size(ind);

% get the 1D nodes and weights
[x1D w1D] = uq_get1DNodesAndWeights(numberOf1DNodes, pc_type);


% multi-dimensional nodes
X = zeros(nquad, ndim);
for i = 1 : nquad
   X(i,:) = x1D(ind(i,:));
end
 
% multi-dimensional weights
W = zeros(nquad, 1);
for i = 1 : nquad

   prodw = 1;
 
   for j = 1 : ndim
        prodw = prodw * w1D(ind(i, j));
   end

   W(i) = prodw;

end

quadrature.nodes   = X;
quadrature.weights = W;
quadrature.nquad = nquad;
