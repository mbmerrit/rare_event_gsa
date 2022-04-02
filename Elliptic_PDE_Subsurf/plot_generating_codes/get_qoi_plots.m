function get_qoi_plots(pde_data, KLa)
nKL = KLa.nKL;
N = 20; 
Xi = randn(N, nKL);
[U xmesh] = get_qoi_evals(Xi, KLa, pde_data);

close all;
plot(xmesh, U, 'linewidth', 2);
xlabel('x');
ylabel('f(x)');
set(gca, 'fontsize', 20);
print('-depsc', 'spe_qoi_rlz');

