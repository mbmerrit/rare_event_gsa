function [qoi] = qoi_wKLE(xi, corrx, corry, sigma)
% This code solve the Darcy flow PDE, the ODE in time, and returns the
% hitting time for a particle with a given starting position. Needs the RVs
% to define the KLE, the KL struct, and data defining the FE mesh for the
% PDE. Implemented by Mike Merritt, March 2021. 
fixpaths

% Get the PDE geometry, mesh, and boundary data
%nx = 80; ny = 40; corrx = 1/2; corry = 1/4;
pde_data = get_pde_data(nx, ny);

% setup the source term in the PDE 
L = 0.005;
% size of point-source at each injection well 
w = [1 1 1 1];
wells = {[0 .2], [0 .4], [0 .6], [0 .8]};
nwells = length(wells);
for i = 1 : nwells
   x0 = wells{i};
   del{i} = @(x,y)( 1/(2*pi*L) * exp(-(1.0/(2*L)*((x-x0(1)).^2 + (y-x0(2)).^2)))); 
end
frhs = @(x, y)(w(1)*del{1}(x, y) + w(2)*del{2}(x, y) + w(3)*del{3}(x, y) + w(4)*del{4}(x, y));
pde_data.frhs = frhs;

% initialize SPE data
XY = pde_data.p';
cd perm_data
%
% specify the vertical slice (SPE data is 3D, we take a 2D slice at a given level) 
level = 32;
init_perm_func(level);
cd('..');
for i = 1 : length(XY)
   mu(i) = perm_func(XY(i,1),XY(i,2));
   %disp(i);
end

% prepare the permeability data
mu = 0.5 * mu(:);   % scale
theta = 1e-4;       % smoothing level (larger ===> smoother)
mu_smooth = get_smooth_mu(mu,pde_data, theta); 

parameters.log_kappa_ave     = mu_smooth;
%parameters.log_kappa_stddev  = .2*max(mu);
approx_mu_inf_norm = 8;
parameters.log_kappa_stddev  = .2*approx_mu_inf_norm;
p = pde_data.p;
e = pde_data.e;
t = pde_data.t;
parameters.corrx = corrx;              % correlation length in x direction
parameters.corry = corry;              % correlation length in y direction
parameters.nx = pde_data.nx; 
parameters.ny = pde_data.ny; 
parameters.mesh = p'; 

nKL = 130; % pick a lot of KL modes
% This function solves the generalized eigenvalue problem for the KLE
KLa_hires = KL_expansion(parameters, nKL);

ntrunc = 126;      % r(k) > 90% 
KLa = truncate_KLa(KLa_hires, ntrunc);    % we use this in computations
%disp('KL expansion set up')

% Now the PDE Solve
[u, ux, uy, c,log_kappa,q] = forward_solve(sigma*xi, KLa, pde_data);

kappa_mid = exp(pdeintrp(pde_data.p,pde_data.t,log_kappa)); % kappa at FE triangle midpoints
v_darcy_t = -[ux;uy] .* kappa_mid;
v_darcy_p = pdeprtni(pde_data.p,pde_data.t,v_darcy_t); % interp velocity at FE nodes

% Solve for several paths through the domain and return hitting times
hitting_times = [];
x0_y = linspace(.3, .7, 1e2);
for i = 1:length(x0_y)
    x0 = [0; x0_y(i)];
    [T, X] = ode_solve(pde_data.p, pde_data.t, x0, v_darcy_p); 
    hitting_times = [hitting_times, T(end)]; % may want to preallocate
end
qoi = hitting_times;

% [T, X] = ode_solve(pde_data.p, pde_data.t, [-1.0, .5], v_darcy_p); 
% qoi = T(end);
end
