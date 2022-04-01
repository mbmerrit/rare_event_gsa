function expf = uq_exp(pcdata, f)
% function: uq_exp.m
% purpose:  computes galerkin exponential of a PCE
%
% usage: expf = uq_exp(pcdata, f)
%
% input:
%    
%    pcdata: struct containing basic PC parameters 
%    f  : is the input PCE
%
% output:
%    expf  : is the exponential of f
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



% The idea here is to solve df/dx=f using forward Euler
% method.

% u is the initial guess of f
x = zeros(size(f));
u = x;
x(1) = f(1);
u(1) = exp(x(1));

% Number of integration steps
nsteps = 2^8;

% Step size
dx = (f - x) / nsteps;

for i = 1 : nsteps
    udx = uq_product(pcdata, u, dx);
    u = u + udx;
end

expf=u;
