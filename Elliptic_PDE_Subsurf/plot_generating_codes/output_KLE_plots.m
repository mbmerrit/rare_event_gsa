function output_KLE_plots(xmesh, U, KLa, KLp)

close all
Ns = length(U);
if Ns < 1000
   U0 = U;
   warning('sample size might be small');
else
   U0 = U(1:1000, :);     % pick the first 1000 for spectrum plot
end

figure(1); 
KLp0 = init_KLE_Pressure(xmesh, 30, U0);
v = KLp0.Dk / KLp0.Dk(1);
semilogy(v, '-o', 'linewidth', 2)
xlabel('k');
ylabel('\lambda_k / \lambda_1');
set(gca, 'fontsize', 20);
print('-depsc', 'spe_output_spectrum');

if Ns >= 5000 
   styl = {'r--o', 'g--o', 'b--o', 'r-o', 'b-o', 'k-o'}; 
   % do a conv study if at least 5000 samples are in U
   nkl = [50 100 500 1000 2500 5000]
   n_nkl = length(nkl);
 
   figure(2); hold on;
   for i = 1 : n_nkl
      Ui = U(1:nkl(i), :);
      KLpi = init_KLE_Pressure(xmesh, 30, Ui);
      v = KLpi.Dk / KLpi.Dk(1);
      semilogy(v, styl{i}, 'linewidth', 2)
   end
   xlabel('k');
   ylabel('\lambda_k / \lambda_1'); 
   set(gca, 'fontsize', 20);
   set(gca, 'yscale', 'log');
   axis tight
   legend('N_{MC} = 50', 'N_{MC} = 100', 'N_{MC} = 500', 'N_{MC} = 1000', 'N_{MC} = 2500', 'N_{MC} = 5000'); 
   print('-depsc', 'spe_output_spectrum_conv');
end


figure; hold on;
u = KLa.sv.^2 ./ KLa.sv(1);
v = KLp.Dk ./ KLp.Dk(1);
semilogy(u, '-o', 'linewidth',2)
semilogy(v, '-o', 'linewidth',2)
xlabel('k');
ylabel('normalized eigenvalues');
set(gca, 'fontsize', 20);
set(gca, 'yscale', 'log');
axis tight
xlim([0 30]);
legend('\lambda_k(C_p)', '\lambda_k(C_q)', 'location', 'best')
grid on;
print('-depsc', 'spe_input_output_spectrum');
