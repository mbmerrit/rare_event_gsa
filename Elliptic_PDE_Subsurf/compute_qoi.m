function [qoi_tau] = compute_qoi(X, method, varargin)

if strcmp(method, 'vary_kle') == 1 % need to re-solve eigenvalue problem each time
    N = size(X,1);
    qoi_tau = zeros(100,N);

    parfor i = 1 : N

        % show progress every set number of samples
        if mod(i, 1e3)==0
            disp(i);
        end
        
        % PDE and ODE evaluation 
        [qoi_tau(:,i)]= qoi_wKLE(X(i,:)', X(i,1), X(i,2), X(i,3));   


    end
elseif strcmp(method,'fix_kle') == 1
    N = size(X,1);
    qoi_tau = zeros(1,N);
    % Initialize with user defined data
    x0 = [0, 0.5]; % position from which the particle is released
    
    if nargin > 2
        KLa = varargin{1}; % If you provide a full KL struct, we can use it
        pde_data = varargin{2}; % also need PDE data
        sigma = varargin{3};
    else
        nx = 50; ny = 50; corrx = 1/2; corry = 1/2;
        % this is the most expensive step! Solve eig problem using Nystrom's method
        [KLa, pde_data] = build_kle_fromparams(nx, ny, corrx, corry);
        sigma = 0.8;
    end
    % improves performance by swithing between series and parallel
    if N == 1
        parforArg = 0;
    else
        parforArg = Inf;
    end
    
    parfor (i = 1 : N, parforArg)   
        % show progress every set number of samples
        if mod(i, 1e4)==0
            disp(i);
        end 
        % Solve the PDE and then the ODE, find the hitting time
        [u, ux, uy, c,log_kappa,q] = forward_solve(sigma*X(i,:)', KLa, pde_data);
        kappa_mid = exp(pdeintrp(pde_data.p,pde_data.t,log_kappa)); % kappa at FE triangle midpoints
        v_darcy_t = -[ux;uy] .* kappa_mid;
        v_darcy_p = pdeprtni(pde_data.p,pde_data.t,v_darcy_t); % interp velocity at FE nodes
        [Time, Xpos] = ode_solve(pde_data.p, pde_data.t, x0, v_darcy_p); 
        qoi_tau(i) = Time(end);
    end
    
%qoi_tau = exp(X); % tail estimation of log-normal example
%mu = [1 2 3 4 5]; sigma = [10 8 6 4 2]; % hyperparameters
%qoi_tau = sum(X .* sigma + mu,2) / sqrt(5); % example from Sehic paper
  
end