function data = get_pde_data(nx, ny)
% Returns PDE geometry and mesh specs in a struct

%load bg_mixed;
load bg_halfdomain_topbottom_neumann;

data.g = g; data.b = b;
if nargin == 0
   nx = 100;
   ny = 50;
   %nx = 80;
   %ny = 40;
end

% Square physical domain.
%[p,e,t]=initmesh(g,'hmax',0.01);
[p,e,t] = poimesh(data.g,nx,ny); % this generates a rectangular mesh
data.p = p;
data.e = e;
data.t = t;
data.nx = nx;
data.ny = ny;

