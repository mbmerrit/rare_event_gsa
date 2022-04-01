function v = perm_func(X,Y)
global XX YY Anomaly logAnomaly

v = interp2(XX,YY,logAnomaly,X,Y,'nearest',0);

