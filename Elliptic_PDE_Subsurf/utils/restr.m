
function [Qu,x] = restr(pde_data, u)


p = pde_data.p;

% defining the QOI as linear functional of solution.
idx      = (p(2,:)==1);      % pick dofs at the boundary y=1

x = p(1, idx);
Qu = u(idx); 
