function [Ti I] = uq_getTotalSensitivity(pcdata, pce, i)
% function: uq_getTotalSensitivity.m
%
% usage: Ti = uq_getTotalSensitivity(pcdata, pce, i)
%
% purpose: returns the total Sobol sensitivity index T_i
% input: 
%   pcdata --- struct containing basic PC info
%   pce    --- vector of PC coefficients
%   i      --- specifies the ith input variable
%              with respect to which we want to 
%              compute global sensitivity.
%              i should be in {1, ..., ndim}
% output:
%   Ti     --- total sensitivity index T_i
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



% extract info from pcdata
ndim = pcdata.ndim;
nord = pcdata.nord;
psiMultiDSq = pcdata.psiMultiDSq;
multiIndex = pcdata.multiIndex;

% get the index set corresponding to \xi_i
I = find(multiIndex(:, i) > 0);

% compute T_i
u = pce.^2;
v = psiMultiDSq;
Ti = dot(u(I), v(I)) / dot(u(2:end), v(2:end));
