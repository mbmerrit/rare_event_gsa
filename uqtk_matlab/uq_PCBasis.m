function Psi = uq_PCBasis(pcdata, xi)
% function: uq_PCBasis.m
% purpose:  return the multidimensional polynomials \Psi_k
%           evaluated at xi.
%
% usage: Psi = uq_PCBasis(pcdata, xi)
% 
% input paremeters:
% 
%    pcdata : struct containing basic PC info
%    xi   : the evaluation point 
%
% output:
%
%    Psi  : a vector with Psi(k) = \Psi_k(xi), k=1, ..., 1+P
%
% note: here the modes of the PC expansion are 
%       numbered 1, ..., P+1 because Matlab array 
%       index must begin with 1. 
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


% check input
if length(xi) ~= pcdata.ndim
   fprintf('error: input vector must be %i dimensional\n', pcdata.ndim);
   Psi = 0; 
   return
end

% get the little psi's: 
psi = uq_psi(xi, pcdata.nord, pcdata.pc_type);

% get the mult-dim Psi's
ndim = pcdata.ndim;
multiIndex = pcdata.multiIndex;
nPCTerms = pcdata.nPCTerms;
Psi = ones(1, nPCTerms);
for i = 1 : ndim
    Psi = Psi .* psi(i, multiIndex(1 : nPCTerms, i) + 1);
end
