function Upc = uq_evalpce(pcdata, pce, xi)
% function: uq_evalpce.m
% purpose: Given the PC coefficients from a (truncated) 
%          PC expansion of random variable, this function 
%          evaluates the PC expansion at the point xi.
%
% usage:   Upc = uq_evalpce(pcdata, pce, xi) 
%
% 
% Parameters:
%    Input: 
%       pcdata --- struct containing basic PC params
%       pce    --- vector of PC coefficients
%       xi     --- point of evaluation
%
%    Output:
%       Upc    --- value of the PC expansion at xi
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

Psi = uq_PCBasis(pcdata, xi);
Upc = dot(pce, Psi);
