function pcdata = uq_pcset(nord, ndim, pc_type)
% function: uq_pcset 
% purpose: sets up all the basic PC parameters
%
% usage: pcdata = uq_pcset(nord, ndim, pc_type)
%
% input: 
% 
%    nord: order of the PCE
%    ndim: the stochastic dimension
%
% output: 
%
%    pcdata: struct containing basic PC information:
%       pcdata.ndim   ---   stochastic dimension
%       pcdata.nord   ---   PCE order
%       pcdata.nclp   ---   number of 1D quadrature points
%       pcdata.pc_type   ---   the choice of PC basis:
%
%                              'LEGENDRE'  ---  Legendre polynomials
%                              'HERMITE'   ---  Hermite polynomials
%                              'LAGUERRE'  ---  Laguerre polynomials
%                              
%       pcdata.x      ---   1D quadrature nodes
%       pcdata.w      ---   1D quadrature weights
%       pcdata.psi    ---   1D orthogonal polynomials at quadrature nodes
%       pcdata.multIndex    ---   the multi-indicies 
%       pcdata.nPCTerms    ---   number of terms in PCE
%       pcdata.psiMultiDSq ---   vector of mean squares of the multi-dim 
%                           polynomials
%
%       sparse storage for triple products 
%          pcdata.nsparse, pcdata.isparse, pcdata.jsparse, pcdata.csparse
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


pcdata.ndim = ndim;
pcdata.nord = nord;
pcdata.nclp = 2 * nord + 1;          % done as in C++ UQToolkit 
pcdata.pc_type = pc_type;

% computing 1D quadrature points and weights
[x, w] = uq_get1DNodesAndWeights(pcdata.nclp, pc_type);
pcdata.x = x;
pcdata.w = w;

% one-dimensional polynomials at collocation points 
pcdata.psi = uq_psi(pcdata.x, nord, pc_type);

% multi-indices and the number of terms in PCE
[pcdata.multiIndex, pcdata.nPCTerms] = uq_initMultiIndex(nord, ndim);

% norms squares
pcdata.psiMultiDSq = uq_evalBasisNormsSquared(nord, ndim, pc_type);

% compute triple moments
pcdata = uq_settripleprod(pcdata);
