
% specify the vertical slice (SPE data is 3D, we take a slice at 
% a given level for 2D test problem), see SPE10 webpage for more info
%level = 50;

% this needs to be run once 
init_perm_func(level);



% the domain in [0,2.2] x [0, 1.2]
x = -1 : 0.01 : 1; 
y = 0 : 0.01 : 1;


% make a plot of log-permeability for illustration
[X Y] = meshgrid(x, y);
K = perm_func(X, Y);
surf(x, y, K); 
shading flat;
view(0,90); 
axis equal tight;


