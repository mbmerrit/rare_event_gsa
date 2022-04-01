function [x, w] = uq_get1DNodesAndWeights(nclp, pc_type)
% function: uq_get1DNodesAndWeights.m
% purpose:  returns the 1D quadrature nodes and weights of the
%           requested Gauss quadrature 
%
% usage: [x, w] = uq_get1DNodesAndWeights(nclp, pc_type)
%
% input:
%
%    nclp     : number of 1D nodes
%    pc_type  : the type of PC basis used; this determines which
%               quadrature rule to use. 
%
% output:
%    x, w     : quadrature nodes and weights  
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


% computing 1D quadrature points and weights
switch pc_type
    case 'LEGENDRE'                      % Gauss-Legendre quadrature
        [x, w] = qrule(nclp, 1, 0, 0);
        w = 0.5 * w;
    case 'HERMITE'                       % Gauss-Hermite quadrature
        [x, w] = qrule(nclp, 9, 0, 0);
        x = x * sqrt(2);
        w = w / sqrt(pi);
    case 'LAGUERRE'                      % Gauss-Laguerre quadrature
        [x, w] = qrule(nclp, 8, 0, 0);
    otherwise
        disp('Unknown quadrature rule');
        return
end

