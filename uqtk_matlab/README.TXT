Introduction:
-------------
This folder contains the source codes for a previous release of the UQ Toolbox.  
The latest version of the UQ Toolkit can be found at https://www.sandia.gov/uqtoolkit/.
UQ Toolbox is a set of Matlab routines for doing UQ and is a part of 
UQ Toolkit. This collection of routines are a result of collaborations
of the following people:

Alen Alexanderian
Maher Salloum
Bert Debusschere
Omar Knio


Basic instructions for the library:
-----------------------------------
The routines for UQ Toolbox are in src_matlab 
folder. To use these routines, either use an 
appropriate addpath statement at the top of each
m-file or put an appropriate addpath statement 
in Matlab's startup file. 


Examples:
---------
Under examples_matlab, there are a few folders containing
some basic illustrations of the routines. 


A brief description of the routines:
------------------------------------
The following are the files for the Matlab UQ Toolbox:

Core routines:
-------------------------------------------------------
uq_pcset.m                ---  sets up basic PC parameters
uq_initMultiIndex.m       ---  sets up multi-indices
uq_get1DNodesAndWeights.m ---  set up 1D nodes and weights
                               (used e.g. for computing triple moments)

PC basis 
-------------------------------------------------------
uq_psi.m                   ---  computes 1D basis polynomials
uq_PCBasis.m               ---  computes mult-D basis polynomials  
uq_evalBasisNormsSquared.m ---  computes 2nd moments of Psi_k 
uq_settripleprod.m         ---  computes triple moments

Operations with PC expansions
------------------------------------------------------
uq_product.m         ---  computes Galerkin products
uq_div.m             ---  computes Galerkin ratios 
uq_exp.m             ---  computes exponential of a PCE
uq_log.m             ---  computes log of a PCE

Note: once we have exp and log we can get 
   * square root : sqrt(pce) = exp(0.5 * log(pce))
   
   * absolute value: abs(pce) = sqrt(pce * pce), where pce*pce is
   computed via uq_product

NISP
------------------------------------------------------
uq_quadtable.m       --- sets up lookup table of quadrature nodes
                         (fully tensorized case for now)
uq_quadrature.m      --- returns mult-dimensional integration nodes and weights
uq_getNISP.m         --- returns the NISP matrix, used to project a vector of
                         data to a PC basis
Sensitivity analysis
------------------------------------------------------
uq_getTotalSensitivity.m  --- computes total Sobol sensitivity indices
uq_getSensitivity.m       --- computes 1st or 2nd order Sobol sensitivity indices
uq_getSensitivityTable    --- computes a matrix containing 1st and 2nd order Sobol indices 

Post-processing
------------------------------------------------------
uq_sample.m          ---  samples a PCE 
                          note: uq_sample uses parfor to speed up the sampling
                          if Matlab's Parallel Processing toolbox is available, 
                          invoke matlabpool before running uq_sample to take advantage
                          of this.
uq_evalpce.m         ---  evaluates a PCE at a point
uq_meanvar.m         ---  computes mean and var of a time series
uq_pdf.m             ---  generates PDF using kde

Additional routines:
------------------------------------------------------
qrule.m              ---  computes 1D quad nodes, weights for Gauss-quadrature
kde.m                ---  kernel density estimation routine
randexp.m            ---  simple routine to get standard exponenetial variables
