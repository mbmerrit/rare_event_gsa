function init_perm_func(lev)

load spe_perm.dat
global XX YY Anomaly logAnomaly;

n1 = 60*85*220;
x = linspace(-1,1,220);
y = linspace(0,1,60);
[XX,YY] = meshgrid(x,y);
A = spe_perm';
A = A(:);
A1 = reshape(A(1:n1),60,220,85);
Anomaly = A1(:,:,lev);
logAnomaly = log(Anomaly);
