function [psiMultiDSq, nPCTerms] = uq_evalBasisNormsSquared(nord, ndim, pc_type)
% function: uq_evalBasisNormsSquared.m
% purpose:  return the 2nd moments of the 
%           multidimensional polynomials \Psi_k.
%
% usage: [psiMultiDSq, nPCTerms] = uq_evalBasisNormsSquared(nord, ndim, pc_type)
%
% input paremeters:
%
%    nord : order of PC expansion
%    ndim : stochastic dimension
%    pc_type : kind of PC basis 
%
% output:
%
%    psiMultiDSq: square norms <Psi_k^2> for k = 0, ..., nPCTerms-1 
%    nPCTerms     : number of terms in the PCE
%
% note: analytical experssions are used for the 2nd moments
%       of 1D polynomials
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

if (nargin < 3)
   pc_type = 'LEGENDRE';   % default
end

[multiIndex, nPCTerms] = uq_initMultiIndex(nord, ndim);
psiMultiDSq = zeros(nPCTerms,1);

switch pc_type

   case 'LEGENDRE' % Legendre polynomials
      norm1d = @(n)( 1.0 / (2*n + 1.0) );
   case 'HERMITE' % Hermite polynomials
      norm1d = @(n)( factorial(n) );
   case 'LAGUERRE' % Laguerre polynomials
      psiMultiDSq = ones(nPCTerms,1);   % norm squares are ones 
      return      
   otherwise
      disp('Unknown quadrature rule');
      return
   end

   % compute norm squares of Psi using the norms squares of 
   % 1D polynomials
   for k = 1 : nPCTerms
      pprod = 1.0;
      for j = 1 : ndim
         pprod = pprod * norm1d(multiIndex(k,j));
      end
      psiMultiDSq(k) = pprod;
   end
end

