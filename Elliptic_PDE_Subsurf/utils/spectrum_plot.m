function spectrum_plot(lambda)

figure(1);
semilogy(lambda./lambda(1), '-o', 'linewidth',2);

figure(2);
ratio = cumsum(lambda) ./ sum(lambda); 
semilogy(ratio, '-o', 'linewidth',2);




