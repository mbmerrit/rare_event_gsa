function ind = uq_quadtable(nclp, ndim)
% function: uq_table
% purpose:  returns a lookup table for integration nodes
%           in a fully tensorized Gauss quadrature grid. 
%
% usage: ind = uq_quadtable(nclp, ndim)
%
% input:
%    nclp: number of samples per dimension
%    ndim: number of stochastic dimensions
%
% output:
%    ind: lookup table for integration nodes 
%            ind(m, i) is the ith-coordinate of the m-th realization.
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


nquad = nclp^ndim;

ind = zeros(nquad, ndim);

for i = 1 : ndim
   ind(1,i) = 1;
end

for m = 2 : nquad

   for i = 1 : ndim
      ind(m,i) = ind(m-1,i);
   end

   nflag = 0;

   ind(m,1) = ind(m,1) + 1;

   if(ind(m,1) > nclp) 
      nflag = 1;
      ind(m,1) = 1;

      for i = 2 : ndim
         if nflag > 0
            nflag = 0;
            ind(m,i) = ind(m,i)+1;

            if ind(m,i) > nclp
               ind(m,i) = 1;
               nflag = 1;
            end

         end
      end
   end
end
