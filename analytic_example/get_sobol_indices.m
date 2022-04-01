function [SobolIndices] = get_sobol_indices(YA_full,YB_full,YC_full)
ndim = size(YC_full, 3);
[N,M] = size(YA_full);
T = zeros(M,ndim);
S = zeros(M,ndim);

if M > 1
for i = 1:M
    % Get vectors of quantities for each M
    YA = YA_full(:,i);
    YB = YB_full(:,i);
    YC = squeeze(YC_full(:,i,:));
    
    % Estimate mean and variance
    muYA = mean(YA);
    muYB = mean(YB);
    Var_Y = (1/2)*mean( (YA-muYA).^2 + (YB-muYB).^2 );

    % Estimate Sobol' indices
    for k = 1:ndim
        %S(i,k) = mean( YB.*(YC(:,k)-YA) )/Var_Y;
        %T(i,k) = (1/2)*mean( (YA - YC(:,k)).^2 )/Var_Y;
        SobolIndices.S_Saltelli(i,k) = mean( YB.*(YC(:,k)-YA) )/Var_Y; % this is our default
        %SobolIndices.S_Sobol(i,k) = (mean( YA.*YD(:,k) ) - mean([muYA muYB])^2 )  /Var_Y;
        SobolIndices.S_Jansen(i,k) = (Var_Y - (1/2)*mean( (YB - YC(:,k)).^2 )) / Var_Y; 
        SobolIndices.T_Saltelli(i,k) = (1/2)*mean( (YA - YC(:,k)).^2 )/Var_Y; % this is our default
        SobolIndices.T_Sobol(i,k) = mean( YA.*( YA - YC(:,k) ) ) /Var_Y;
        SobolIndices.T_Homma(i,k) = (Var_Y - mean(YA.*YC(:,k)) + mean([muYA muYB])^2 ) / Var_Y;
    end
end  
else
    YA = YA_full;
    YB = YB_full;
    YC = YC_full;
    muYA = mean(YA_full);
    muYB = mean(YB_full);
    Var_Y = (1/2)*mean( (YA-muYA).^2 + (YB-muYB).^2 );

    % Estimate Sobol' indices COMPARING ESTIMATORS
    for k = 1:ndim
        %SobolIndices.S(k) = mean( YB.*(YC(:,k)-YA) )/Var_Y;
        %SobolIndices.T(k) = (1/2)*mean( (YA - YC(:,k)).^2 )/Var_Y;
        SobolIndices.S_Saltelli(k) = mean( YB.*(YC(:,k)-YA) )/Var_Y; % this is our default
        %SobolIndices.S_Sobol(k) = (mean( YA.*YD(:,k) ) - mean([muYA muYB])^2 )  /Var_Y;
        SobolIndices.S_Jansen(k) = (Var_Y - (1/2)*mean( (YB - YC(:,k)).^2 )) / Var_Y; 
        SobolIndices.T_Saltelli(k) = (1/2)*mean( (YA - YC(:,k)).^2 )/Var_Y; % this is our default
        SobolIndices.T_Sobol(k) = mean( YA.*( YA - YC(:,k) ) ) /Var_Y;
        SobolIndices.T_Homma(k) = (Var_Y - mean(YA.*YC(:,k)) + mean([muYA muYB])^2 ) / Var_Y; 
    end
end

        
end