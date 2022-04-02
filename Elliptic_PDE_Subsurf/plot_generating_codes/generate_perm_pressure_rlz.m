% make sure to run fixpaths in the main folder for the
% codes for this problem. Run the code from the folder
%  gsa_asm/codes/Elliptic_PDE_Subsurf
%load plotting_data
clear u log_kappa;
NKLa = length(KLa.sv);
rng(2345);
Ns = 4;
X = randn(Ns, NKLa); 
close all
for i = 1 : Ns
   [u{i} log_kappa{i}] = forward_solve(X(i, :), KLa, pde_data);

   mu(i)  = min(u{i}(:));
   Mu(i)  = max(u{i}(:));
   
   ma(i)  = min(log_kappa{i}(:));
   Ma(i)  = max(log_kappa{i}(:));

end

mu = min(mu); Mu = max(Mu);
ma = min(ma); Ma = max(Ma);
cbar_flag = 0;
img_w = 24;
img_h = 12;
for i = 1 : Ns
   figure(1+i); 
   plot_vec(pde_data, u{i}, 3, cbar_flag) 
   caxis([mu Mu]);
   axis off
   if i == Ns
      orig_size = get(gca, 'Position');
      colorbar;
      set(gca, 'Position', orig_size);
   end
   fname = ['spe_pressure_' num2str(i)];
   figsize(img_w,img_h,'cm')
   print('-dpng', fname);
 
   figure(Ns+i); 

   plot_vec(pde_data, log_kappa{i}, 1, cbar_flag); view(2)
   caxis([ma Ma]);
   axis off
   if i == Ns
      orig_size = get(gca, 'Position');
      colorbar;
      set(gca, 'Position', orig_size);
   end
   zlim([ma Ma]); 
   fname = ['spe_perm_' num2str(i)];
   figsize(img_w,img_h,'cm')
   print('-dpng', fname);
   
end
