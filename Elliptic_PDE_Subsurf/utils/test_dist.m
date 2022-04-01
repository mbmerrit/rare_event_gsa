function test_dist(U, U2, idx)
[f xi] = ksdensity(U(:,idx), 'support', 'positive', 'npoints', 1e3);
[f2 xi2] = ksdensity(U2(:,idx), 'support', 'positive', 'npoints', 1e3);
plot(xi, f, 'linewidth',2)
hold on;
plot(xi2, f2, 'linewidth',2)
legend('full', 'reduced');
