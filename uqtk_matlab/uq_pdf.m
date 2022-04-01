function [xmesh, density] = uq_pdf(data, nmesh, lb, ub)
% function: uq_pdf.m
% purpose:  uses kde to approximate the pdf of a r.v. 
%           sampled in vector data.
%
% Usage: [xmesh, density] = uq_pdf(data, nmesh, lb, ub)
%
% input:
%    data    --- sampling of the r.v.
%    nmesh    --- number of points in the interval for
%                 PDF approximation.
%    lb, ub (optional)
%    if lb and ub are absent, uses min and max 
%    of data as the pdf interval, if ub is absent
%    makes the pdf over:
%
%        [mean(data) - lb * std(data), mean(data) + lb*std(data)]
%
%    if both lb and ub are specified:
%
%        [mean(data) - lb*std(data), mean(data) + ub*std(data)]
% 
% output: 
%    xmesh, density
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



avg = mean(data);
sigma = std(data);

if (nargin < 3)
   m = min(data);
   M = max(data);
elseif (nargin < 4) 
   m = avg - lb*sigma;
   M = avg + lb*sigma;
else
   m = avg - lb*sigma;
   M = avg + ub*sigma;
end

% KDE estimation of the PDF:
[bandwidth, density, xmesh] = kde(data, nmesh, m, M);
