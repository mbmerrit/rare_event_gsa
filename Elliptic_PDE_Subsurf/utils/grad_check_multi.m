function [err] = grad_check_multi(G, f, xi, KLa, KLp, pde_data, k)
% k is the kth gradient column/eval, this function is very similar to it's
% predesecor grad_check, except this outputs the error, and it excepts a
% vector of indices k
 
ndim = length(KLa.sv);
f = f(k);
G = G(:,k);

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



