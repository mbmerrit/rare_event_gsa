function [S I] = uq_getSensitivity(pcdata, pce, i, j)
% function: uq_getSensitivity.m
%
% usage: S = uq_getSensitivity(pcdata, pce, i, j)
%
% purpose: returns the first or second order sensitivity index.
%
% input: 
%   pcdata --- struct containing basic PC info
%   pce    --- vector of PC coefficients
%   i,j    --- specifies the ith and jth input variable
%              with respect to which we want to 
%              compute global sensitivity.
%              i,j should be in {1, ..., ndim}
%          
%          If only one index i is provided, the returned value will be
%          the first order sensitivity index S_i, whereas if both i and
%          j are supplied, the return value will be the second order sensitivity
%          index S_ij.
% output:
%   S     --- first order index S_i if only i is provided, 
%             or second order index S_ij if i and j are supplied.
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

% extract info from pcdata
ndim = pcdata.ndim;
nord = pcdata.nord;
psiMultiDSq = pcdata.psiMultiDSq;
multiIndex = pcdata.multiIndex;


% get the appropriate index set for the sensitivity index:
idx = 1 : ndim;
if nargin < 4          % only i is provided
   excluded_dims = multiIndex(:, idx ~= i);
   I = find( (multiIndex(:,i) > 0) & (sum(excluded_dims, 2) == 0) )
else                   % both i, j are provided
   excluded_dims = multiIndex(:, idx ~= i & idx ~= j);
   I = find( multiIndex(:, i) > 0 & multiIndex(:, j) > 0 & sum(excluded_dims, 2) == 0);
end

% compute the sensitivity index
u = pce.^2;
v = psiMultiDSq;
S = dot(u(I), v(I)) / dot(u(2:end), v(2:end));
