function pce = uq_div(pcdata, pce1, pce2)
% function: uq_div.m
% purpose:  computes galerkin division of two PCEs
%
% usage: pce = uq_div(pcdata, pce1, pce2)
%
% input:
%
%    pcdata: struct containing basic PC parameters
%    pce1, pce2: these are the pc modes of the expansion
%                we are dividing
%
% output:
%    pce: will contain the  pce1/pce2
% 
% basic idea: Given PCE of random variables u, v, 
%             we want to find w such that 
%             w = u / v
%             or alternatively, w * v = u.
%             A weak formulation of the above 
%             and using Galerkin products allows us
%             to formulate the problem as 
%               
%                 A [w] = [u]                 
%            
%             where [u] denotes a vector containing
%             the PC modes of u (stored in pce1).
%             the entries of the matrix A are in terms of
%             the triple moments <Psi_i Psi_j Psi_k> and 
%             PC modes of v (stored in pce2).
%             
% note: as usual, sparse storage for triple moments is used
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




%------------------------------------------
% the basic setup
%------------------------------------------
nsparse = pcdata.nsparse;
isparse = pcdata.isparse;
jsparse = pcdata.jsparse;
csparse = pcdata.csparse;
nPCTerms = pcdata.nPCTerms;
psiMultiDSq = pcdata.psiMultiDSq;

%------------------------------------------
% setup the matrix A for the linear system:
%   
%    A * pce = pce1
%------------------------------------------
A = zeros(nPCTerms);
for k = 1 : nPCTerms
   for s = 1 : nsparse(k)
      m1 = isparse(k, s);
      m2 = jsparse(k, s);
      c  = csparse(k, s) / psiMultiDSq(k);
      A(k,m1) = A(k,m1) + pce2(m2) * c;
   end
end

%------------------------------------------
% Solve the system with Gauss elimination
%------------------------------------------
pce = A \ pce1(:);
pce = pce';        % make row vector
