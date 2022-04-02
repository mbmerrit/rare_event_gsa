function param_KLE_plots(pde_data, KLa)

% plot mean log-permeability
figure;
generate_mean_field;
print('-dpng', 'param_KLE_mean');

% plot normalized eigenvalues
figure;
v = KLa.sv.^2./KLa.sv(1)^2;
semilogy(v, '-o', 'linewidth', 2);
axis tight;
xlim([0 1000]);
set(gca, 'fontsize', 20);
xlabel('k');
ylabel('\lambda_k/\lambda_1');
print('-depsc', 'param_KLE_plot');


% plot ratios
figure
hold on;
r = KLa.r;
semilogy(r, '-', 'linewidth', 2);
hold on;
rp = [.25 .5 .75 .9 .99 1];
ridx = [1 5 10 50 126 400 1000];
loglog(ridx, r(ridx), 'ko', 'markerfacecolor', 'k');
set(gca, 'ytick', rp);
set(gca, 'xtick', ridx);
grid on
axis tight;
xlim([0 1000]);
set(gca, 'fontsize', 20);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
xlabel('k');
ylabel('r_k');
print('-depsc', 'param_KLE_ratio_plot');


