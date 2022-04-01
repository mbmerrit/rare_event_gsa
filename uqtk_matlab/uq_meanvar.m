function [mu sigma2] = uq_meanvar(pcdata, pce)
% function: uq_meanvar.m
% purpose:  given a pc expansion (or a collection of them)
%           compute mean and variance. This is useful in 
%           getting time evolution of mean and variance 
%           given a time series of PC coefficients
%           stored as rows of the pce matrix. 
%
% usage: [mu sigma2] = uq_meanvar(pcdata, pce)
%
% input:
%           pcdata          struct containing basic PC info
%         
%           pce             M x nPCTerms matrix, whose each row
%                           is a pc expansion (coefficients)
% 
% outout:
%           mu, sigma2  --- mean and variance, both M x 1
%                           vectors 
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



% setup PC info
psiMultiDSq = pcdata.psiMultiDSq;

% get the mean:
mu = pce(:,1);

% compute variance:
A = pce(:,2:end).^2;
b = psiMultiDSq(2:end);
sigma2 = A * b;

