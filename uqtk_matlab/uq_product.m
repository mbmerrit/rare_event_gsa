function pce = uq_product(pcdata, pce1, pce2)
% function: uq_product.m
% purpose:  computes galerkin product of two PCEs
%
% usage: pce = uq_product(pcdata, pce1, pce2)
%
% input:
%    
%    pcdata: struct containing basic PC parameters 
%    pce1, pce2: these are the pc modes of the expansion
%                we are multiplying
%
% output:
%    pce: will contain the product of pce1 and pce2
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



% the basic setup
nsparse = pcdata.nsparse;
isparse = pcdata.isparse;
jsparse = pcdata.jsparse;
csparse = pcdata.csparse;
psiMultiDSq = pcdata.psiMultiDSq;

% the output vector
pce = zeros(size(pce1));

% perform (sparse) galerkin product
for k = 1 : pcdata.nPCTerms
   sum = 0.0;

   for l = 1 : nsparse(k)

      sum = sum + pce1(isparse(k,l)) * pce2(jsparse(k,l)) * csparse(k,l);

   end

   pce(k) = sum / psiMultiDSq(k);

end
