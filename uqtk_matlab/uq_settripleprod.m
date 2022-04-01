function  pcdata_out = uq_settripleprod(pcdata_in)
% function: uq_settripleprod.m
% purpose: computes the triple products 
%             <Psi_i Psi_j Psi_k> 
%
% usage: pcdata_out = uq_settripleprod(pcdata_in) 
%
% input:
%    pcdata_in: contains basic pc data
%
% output: 
%    pcdata_out: will be a copy of the
%                pcdata_in except the fields
%                for triple products are filled out
%
% note: sparse storage is used for efficiency
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

% these are extracted for convenience
nord = pcdata_in.nord;
nclp = pcdata_in.nclp;
ndim = pcdata_in.ndim;
multiIndex  = pcdata_in.multiIndex;
nPCTerms  = pcdata_in.nPCTerms;
w    = pcdata_in.w;
psi  = pcdata_in.psi;

apow=zeros(nord + 1, nord + 1, nord + 1);

for k = 1 : nord + 1
    for j = 1 : nord + 1
        for i = 1 : nord + 1
            sum = 0;
            for m = 1 : nclp
                sum = sum + psi(m,i)*psi(m,j)*psi(m,k)*w(m);
            end
            apow(i,j,k) = sum;
        end
    end
end



nsparse(1:nPCTerms)=0;

% tolerance to prune out zero coeffs
tol = 1e-4;

for k = 1 : nPCTerms
    for j = 1 : nPCTerms
        for i = 1 : nPCTerms
            aprod = 1;
            for m = 1 : ndim
                l1 = multiIndex(k, m) + 1;
                l2 = multiIndex(j, m) + 1;
                l3 = multiIndex(i, m) + 1;
                aprod=aprod*apow(l1, l2, l3);
            end    
            if (aprod > tol)
                nsparse(k) = nsparse(k) + 1;
                isparse(k, nsparse(k)) = i;
                jsparse(k, nsparse(k)) = j;
                csparse(k, nsparse(k)) = aprod;
            end
        end
    end
end

% setup the output
pcdata_out = pcdata_in;
pcdata_out.nsparse = nsparse;
pcdata_out.isparse = isparse;
pcdata_out.jsparse = jsparse;
pcdata_out.csparse = csparse;
