%  This script is meant to demonstrate how the QoI is computed in the
%  following steps: 1) The KL expansion is set up and the random
%  permeability field is generated. 2) The Darcy flow PDE is solved using
%  the finite element method. 3) The time-evolution ODE is solved for a
%  single particle in the domain. Each of these steps must take place for
%  every QoI evaluation, with many evaluations required for rare event
%  estimation. Implemented March 2021, Mike Merritt. 

%
fixpaths
% Get the PDE geometry, mesh, and boundary data
nx = 50; ny = 50; corrx = .5; corry = .5;
pde_data = get_pde_data(nx, ny);
%
% setup the source term in the PDE 
%
L = 0.005;
% size of point-source at each injection well, this part is not used later 
w = [1 1 1 1]; wells = {[0 .2], [0 .4], [0 .6], [0 .8]}; nwells = length(wells);
for i = 1 : nwells
   x0 = wells{i}; del{i} = @(x,y)( 1/(2*pi*L) * exp(-(1.0/(2*L)*((x-x0(1)).^2 + (y-x0(2)).^2)))); 
end
frhs = @(x, y)(w(1)*del{1}(x, y) + w(2)*del{2}(x, y) + w(3)*del{3}(x, y) + w(4)*del{4}(x, y));
pde_data.frhs = frhs;

%
% initialize SPE data
%
XY = pde_data.p';
cd perm_data
%
% specify the vertical slice (SPE data is 3D, we take a slice at
% a given level for 2D test problem), see SPE10 webpage for more info
% https://www.spe.org/web/csp/datasets/set02.htm
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

nKL = 1000; % pick a lot of KL modes
% This function solves the generalized eigenvalue problem for the KLE
KLa_hires = KL_expansion(parameters, nKL);
disp(['Capturing 90% of the KL variance with ',num2str(find(KLa_hires.r > .9, 1)),' terms'])

ntrunc = 126;      % r(k) > 90% is the goal for effective KL use
KLa = truncate_KLa(KLa_hires, ntrunc);    % we use this in computations
disp('KL expansion set up')

%%  This block does the PDE solve, computes the Darcy velocity field at 
% discrete points and then solves the ODE for particle motion

% Plot PDE properties and solve PDE 
pdemesh(pde_data.p,pde_data.e, pde_data.t); title('finite element mesh')% plot the FE mesh
xi = 1*randn(ntrunc,1); % these are the RVs that characterize the KLE
[u, ux, uy, c,log_kappa,q] = forward_solve(xi, KLa, pde_data);

figure()
y = evaluate_KL_expansion(KLa, xi  );
pdeplot(p,e,t,'XYData',log_kappa, 'colormap','jet'); title('log kappa'); %pbaspect([1 1 1]); 

% create model and populate items
model = createpde(); geometryFromEdges(model,pde_data.g);
% plot the PDE solution, pressure field, then gradients
figure()
% pdeplot(p,e,t,'XYData',u); hold on
pdeplot(p,e,t,'flowdata',-[ux; uy],'XYData',u,'ColorMap','default',...
    'contour','on','flowstyle','arrow'); title('pressure solution, pressure gradient field')
axis tight; %pbaspect([1 1 1]);

figure()
kappa_mid = exp(pdeintrp(p,t,log_kappa)); % kappa at FE triangle midpoints
v_darcy_t = -[ux;uy] .* kappa_mid;
v_darcy_p = pdeprtni(p,t,v_darcy_t); % interp velocity at FE nodes
pdeplot(p,e,t,'flowdata',v_darcy_t,'XYData',u,'ColorMap','default',...
    'contour','off','flowstyle','arrow'); title('Pressure solution, Darcy velocity field')
% Get the X,Y coordinates of the FE cells, there are 2*nx*ny centroids
[XY_centroids] = get_centroid_coordinates(nx, ny, 0, 1, 0, 1);
%%
% Plot several paths through the medium, using the ODE solve
hitting_times = [];
for i = linspace(.01, .99, 50) % plot streamlines from bottom to top
    x0 = [0; i]; hold on
    [T, X] = ode_solve(p, t, x0, v_darcy_p); plot(X(:,1), X(:,2), 'm-')
    hitting_times = [hitting_times, T(end)]; % may want to preallocate
end
[T, X] = ode_solve(p, t, [0,.5], v_darcy_p); plot(X(:,1), X(:,2), 'k-','linewidth',2)

