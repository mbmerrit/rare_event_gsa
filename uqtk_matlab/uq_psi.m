function psi = uq_psi_1D(x, nord, pc_type)
% function: uq_psi
% purpose:  Generates 1-D polynomials
%
% usage: psi = uq_psi(x, nord, pc_type)
%           
% input:
%    x   : column vector of points
%    nord: expansion order
%    pc_type: polynomial pc_type
%
% output:
%    psi:  polynomial evaluated at points x [length(x) * nord]
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


nclp = length(x);
psi = zeros(nclp, nord + 1);

switch pc_type
    case 'LEGENDRE'
        psi(:, 1) = 1;
        
        if (nord > 0)
            psi(:, 2) = x;

            % computing 1D legendre polynomials
            for i = 3 : nord+1
                psi(:,i)= (2*i-3)/(i-1) * x(:) .* psi(:, i-1) - (i-2)/(i-1) * psi(:, i-2);
            end
        end
    case 'HERMITE'
        psi(:, 1) = 1;
        
        if (nord > 0)
            psi(:, 2) = x;

            % computing 1D hermite polynomials
            for i = 3 : nord+1
                psi(:, i) = x(:) .* psi(:, i-1) - (i-2) * psi(:, i-2);
            end
        end
    case 'LAGUERRE'
        psi(:, 1) = 1;
        
        if (nord > 0)
            psi(:, 2) = 1 - x;

            % computing 1D laguerre polynomials
            for i = 3 : nord+1
                psi(:,i)= 1/(i-1) * ( (2*i - 3 - x(:)) .* psi(:, i-1) - (i-2) * psi(:, i-2));
            end
        end
    otherwise
        disp('Error: Unknown polynomial pc_type')
end
        
