function [KLp] = init_KLE_Pressure(xmesh, NKl, U)

%
% prepare mesh data 
%
mesh.nodes = xmesh;  
dx = xmesh(2)-xmesh(1);
% trapz weights
w = ones(size(xmesh))*2; 
w(1) = 1;
w(end) = 1; 
w = 0.5 * dx * w;    
mesh.weights = w(:);
mesh.nx = length(xmesh);

[C,R,Ubar,Uc] = get_pressure_covariance(U);
[Dk Vk] = get_KL_modes(mesh, NKl, C);

% compute the mode samples
Ns = size(U, 1)
for i = 1 : NKl
   sigma_i = sqrt(Dk(i));
   for j = 1 : Ns
      pk(i,j) = (1/sigma_i) * dot(w, (Uc(j,:)' .* Vk(:,i)));
   end
end    
KLp.Dk = Dk;
KLp.Vk = Vk;
KLp.pbar = Ubar;
KLp.pk = pk;

%
%  --------------- Sub-functions ----------------

function [D V] = get_KL_modes(mesh, N_KL_p, C)

Msqrt = diag(sqrt(mesh.weights));
CM = Msqrt * (C * Msqrt);
[V,lambda] = eigs(CM, N_KL_p);
D = diag(lambda);
V = Msqrt \ V;

[D idx] = sort(D, 'descend');
V = V(:, idx);





