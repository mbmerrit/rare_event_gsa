function [u, ux, uy, c,log_kappa,q] = forward_solve(xi, KLa, pde_data)
%
% Testing forward solve:
%    -div(kappa grad u) + alpha * u = q
%
    
p = pde_data.p;
e = pde_data.e;
t = pde_data.t;
b = pde_data.b;
    
Q = pde_data.frhs(p(1, :)', p(2, :)');
q = 0*pdeintrp(p,t, Q(:)); % use zero source
log_kappa = evaluate_KL_expansion(KLa, xi); 
%log_kappa = p(2, :).^2 + 1; 
%log_kappa = log_kappa(:);
kappa = exp(log_kappa);
c = pdeintrp(p,t, kappa); 
alpha = 0;        
[K,F] = assempde(b,p,e,t,c,alpha,q);
u = K \ F;

[ux,uy] = pdegrad(p,t,u);

