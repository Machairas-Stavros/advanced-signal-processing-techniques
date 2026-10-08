function nrmse = normalizedRootMeanSquareError(x_estim,x)
    %This function calculates the normalized root mean square error (NRMSE)
    %between an estimated signal and the original signal
    
    rmse = sqrt(sum((x_estim-x).^2)/length(x));
    nrmse = rmse/(max(x(:))-min(x(:)));
end