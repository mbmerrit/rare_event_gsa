function U = uq_sample(pcdata, pce, Nsamp)
% function: uq_sample.m
% purpose: Given the PC coefficients from a (truncated) PC expansion of 
% random variable, this function samples the PC expansion. 
% 
% usage: U = uq_sample(pcdata, pce, Nsamp)
%
% Parameters:
%    Input: 
%       pcdata: struct containing basic PC parameters
%       pce:   the vector of PC coefficients
%       Nsamp: number of realizations to be used in sampling 
%
%    Output:
%       U:     the PC expansion computed at sample realizations
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


U = zeros(1,Nsamp);

% pick random number generator based on pc_type:

switch pcdata.pc_type
    case 'LEGENDRE'                           
       randsampler = inline('2*rand(1,x)-1');      % Uniform(-1,1)             
    case 'HERMITE'                           
       randsampler = inline('randn(1,x)');         % standard normal
    case 'LAGUERRE'                           
       randsampler = inline('randexp(1,x)');       % Standard Exponential
    otherwise
        disp('Unknown quadrature rule');
        return
end

ndim = pcdata.ndim;


% Note: in a multicore machine can use parfor to speed up this loop
for s = 1 : Nsamp
    % get a random sample point
    xi = randsampler(ndim);

    % get the mult-dim Psi's
    Psi = uq_PCBasis(pcdata, xi);
 
    % compute PC expansion:
    U(s) = dot(pce, Psi);
end
