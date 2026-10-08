function third_order_cumulants = ThirdOrderCumulantsIndirect(sample,K,M,maxshifts,pltspr)
    %This function calculates the 3rd order cumulants for a time series 
    %'sample' using the indirect method given the number of maximum lags 
    %'maxshifts' that are going to be used the parameters K and M for the 
    %spliting of the original sample to subsamples. The 'pltspr' variable 
    %is a logical variable which enables the display of the plot for the
    %surface of the third order cumulants if true. If false or not defined 
    %(automatically takes the value false) and disables the display of the 
    %pre-mentioned plot.
    
    if nargin == 4
        pltspr = false;
    end
    
    % Step 1
    subsamples = NaN([M K]);
    for i=0:K-1
        subsamples(:,i+1) = sample(i*M+1:(i+1)*M);
    end

    % Step 2
    subsamples = subsamples-mean(subsamples);

    % Step 3
    shiftings = -maxshifts:maxshifts;
    r3x_taf1_taf2 = NaN([length(shiftings) length(shiftings) K]);
    for k=1:K
        for t1=1:length(shiftings)
            for t2=1:length(shiftings)
                r3x_taf1_taf2(t1,t2,k) = (1/M)*sum(subsamples(:,k).*arrayShift(subsamples(:,k),shiftings(t1)).*arrayShift(subsamples(:,k),shiftings(t2)));
            end
        end
    end

    % Step 4
    c3x_taf1_taf2 = sum(r3x_taf1_taf2,3)/K;

    switch pltspr
        case true
            figure()
            surf(shiftings,shiftings,c3x_taf1_taf2);
            title('3^{rd} Order Cumulants','Indirect Method: K='+string(K)+', M='+string(M)+', L_3='+string(maxshifts));
            xlabel('\tau_1');
            ylabel('\tau_2');
            zlabel('c^x_3(\tau_1,\tau_2)');
            colorbar;
    end
    
    third_order_cumulants = c3x_taf1_taf2;
end