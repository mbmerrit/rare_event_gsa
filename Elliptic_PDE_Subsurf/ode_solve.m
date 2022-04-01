function [T, X] = ode_solve(p, t, x0, v_darcy)
%
% Given a Darcy velocity field, evaluated at node points p, this code
% solves the ODE system dx/dt = v_darcy

if nargin < 3    % default rate constants
    error('Need an initial condition for the ODE!')
elseif size(v_darcy') ~= size(p)
    error('Need Darcy velocity evaluations at the nodes, dimension error')
end

% For a speedup, you could also build the interpolants outside the solver
F1 = pdeInterpolant(p,t,v_darcy(:,1));
F2 = pdeInterpolant(p,t,v_darcy(:,2));

options = odeset('AbsTol',1e-8,'RelTol',1e-8, 'Events', @hit_boundary);
frhs = @(T, X)(darcy_interp(p, t, X, v_darcy, F1, F2));
[T,X] = ode45(frhs, [0, 1e6], x0, options); % T_final should be inf

%-------------- sub-functions ----------------
function [v_at_x] = darcy_interp(p, t, X, v_darcy, F1, F2)
% Given a generic point in the domain, this computes the Darcy velocity at
% that point given only sparse data at the FE midpoints. This code requires
% the Darcy velocity to be known at the FE nodes, not the midpoints!
% For now, we pass the PDE interpolant. If we not, we need p, t, and v

% The interpolated Darcy velocity is the time derivative of position
v_at_x = [evaluate(F1,X(1),X(2)); evaluate(F2,X(1),X(2))];
if double(isnan(v_at_x)) == ones(2,1)
    v_at_x = zeros(2,1);
end

% Event function, terminates the ODE solver if the particle exits domain
function [value, isterminal, direction] = hit_boundary(T, X)
value      = double(any(X(1) >= 1.0));
isterminal = 1;   % Stop the integration
direction  = 0;


