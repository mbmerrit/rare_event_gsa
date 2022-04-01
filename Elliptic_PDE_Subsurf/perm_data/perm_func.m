function v = perm_func(X,Y)
global XX YY Anomaly

v = interp2(XX,YY,Anomaly,X,Y,'nearest',0);

