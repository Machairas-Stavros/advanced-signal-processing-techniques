function [c3x_t1_t2,impulse_response,signal_x_estim,nrmse] = cumulantsImpulseNRMSE(signal_x,signal_v,q,K,M,L3)
    %This function calculate the 3rd Order Cumulants of a given signal_x
    %estimates the impulse response given the third order cumulants of the
    %signal_x and the non gaussian noise signal_v of the system and
    %calculates the NRMSE between the signal_x and its estimation.
    
    c3x_t1_t2 = ThirdOrderCumulantsIndirect(signal_x,K,M,L3);

    impulse_response = GiannakisFormula(c3x_t1_t2,q);

    signal_x_estim = conv(signal_v,impulse_response);
    signal_x_estim = signal_x_estim(1:length(signal_x));
        
    nrmse = normalizedRootMeanSquareError(signal_x_estim,signal_x);
end