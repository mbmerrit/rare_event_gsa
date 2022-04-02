% make sure to run fixpaths in the main folder for the
% codes for this problem. Run the code from the folder
%  gsa_asm/codes/Elliptic_PDE_Subsurf
%load plotting_data
clear u log_kappa;
NKLa = length(KLa.sv);
Ns = 4;
X = randn(Ns, NKLa); 
close all
cbar_flag = 1;
plot_vec(pde_data, KLa.mean, 2, cbar_flag); view(2);
hold on;
wells = {[-.6 .2], [-.2 .4], [.2 .6], [.6 .8]};
for i = 1 : 4
   x0 = wells{i};
   plot(x0(1), x0(2), 'o', 'markerfacecolor', 'k', 'markersize', 14);
end
alpha .6
axis off
cbar_flag = 0;
img_w = 24;
img_h = 12;
fname = ['spe_perm_mean'];
figsize(img_w,img_h,'cm')
print('-dpng', fname);


