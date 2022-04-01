function plot_vec(pde_data, u, plot_opt, cbar_flag)
%
% purpose: plot a function with PDE-toolbox tools
%    
% input:   pde_data   (struct containing PDE data)
%          u          (function to be plotted)
%          plot_opt   (1=surf, 2=contour without lines, 3=contour with lines)
%          cbar_flag  (0: no colorbar; 1: with colorbar)
%  


p = pde_data.p;
e = pde_data.e;
t = pde_data.t;
b = pde_data.b;

if nargin < 4
   cbar_flag = 1;
elseif nargin < 3
   plot_opt = 1;
   cbar_flag = 1;
elseif nargin < 2 
   error('too few arguments ...');
end

cbar = {'off', 'on'};

if plot_opt == 1

   pdeplot(p,[],t,'XYData',u,'XYStyle','interp',...
      'ZData',u,'ZStyle','continuous','Contour','off',...
      'ColorBar', cbar{cbar_flag+1});
   view(30,20);
elseif plot_opt == 2
   pdeplot(p,[],t,'XYData',u,'XYStyle','interp',...
      'Contour','off',...
      'ColorBar', cbar{cbar_flag+1});
elseif plot_opt == 3
   pdeplot(p,[],t,'XYData',u,'XYStyle','interp',...
      'Contour','on',...
      'ColorBar', cbar{cbar_flag+1});
end 
colormap jet
axis equal tight;

set(gca, 'fontsize', 20);
