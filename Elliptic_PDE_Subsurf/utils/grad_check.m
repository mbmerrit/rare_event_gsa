function grad_check(KLa, KLp, pde_data, k, j)

ndim = length(KLa.sv);
xi = randn(ndim, 1);
[f,G] = get_kl_mode_sensitivity(xi',KLa,KLp,pde_data);

f = f(k);
G = G(:,k);

if nargin < 5
   d = rand(ndim, 1); d = d / norm(d);
else
   d = zeros(ndim, 1); d(j) = 1;
end


% compute directional derivative
d = rand(ndim, 1); d = d / norm(d);

t = 1/2;
for i = 1 : 9
   f1 =  get_kl_mode_sensitivity((xi+t*d)', KLa, KLp,pde_data);
   f1 = f1(k);
   FD_Grad = (f1 - f) / t;
   Adj_Grad = dot(G, d);
   err(i) = abs(FD_Grad - Adj_Grad) / abs(Adj_Grad);
   fprintf('%10.4e\n', err(i));
   t = t / 2;
end
idx = [1 : 9];
semilogy(idx, err, '-o', 'linewidth',2);

