function [f,G] = get_qoi(xi, KLa, KLp, pde_data, k)


% Forward problem:
%    -div(kappa grad u) = q
%

p = pde_data.p;
e = pde_data.e;
t = pde_data.t;
b = pde_data.b;

m = length(xi);
ndof     = size(p,1);

% defining the QOI as linear functional of solution.
idx      = (p(2,:)==1);      % pick dofs at the boundary y=1

xx = p(1, idx);
psi = zeros(ndof, 1);
pbar = zeros(ndof,1);
pbar(idx) = KLp.pbar;

xx = p(1,:)'; yy = p(2,:)';
xmesh = xx(idx);

%L = 0.01;
%frhs{1} = @(x,y)( 2/(2*pi*L) * exp(-(1.0/(2*L)*((x+.6).^2 + (y-.5).^2)))); % point source
%frhs{2} = @(x,y)( 5/(2*pi*L) * exp(-(1.0/(2*L)*((x-0).^2 + (y-.5).^2)))); % point source
%frhs{3} = @(x,y)( 2/(2*pi*L) * exp(-(1.0/(2*L)*((x-.6).^2 + (y-.5).^2)))); % point source
%Q = 0;
%for i = 1 : 3
%   Q = Q + frhs{i}(xx, yy);
%end
Q = pde_data.frhs(xx, yy);
q = pdeintrp(p,t, Q(:));

% forward
log_kappa = evaluate_KL_expansion(KLa, xi);
kappa = exp(log_kappa);
c=pdeintrp(p,t, kappa); 
alpha = 0;        

[K,F]=assempde(b,p,e,t,c,alpha,q);

% forward solve 
u = K\F;

% compute QoI
fac  = 1/sqrt(KLp.Dk(k));
psi(idx) = KLp.Vk(:, k);
% compute QoI
f = fac * trapz(xmesh, (u(idx)-pbar(idx)).*psi(idx));

W = diag(get_trapz(xmesh));
v = W * psi(idx);
qadj = zeros(size(psi));
qadj(idx) = v;

if nargout == 2
   % adjoint solve 

   phi = K \ -(fac*qadj);

   % sensitivities
   E = KLa.basis;
   sv = KLa.sv;
   for j=1 : m
      coeff = pdeintrp(p,t, sv(j) * E(:,j) .* kappa);
      [K,~] = assempde(b,p,e,t,coeff,alpha,1);
      G(j)  = phi'*(K*u);
   end
   G = G(:);
end
