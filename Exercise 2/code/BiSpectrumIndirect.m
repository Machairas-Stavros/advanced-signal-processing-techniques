function bispec = BiSpectrumIndirect(sample,K,M,maxshifts,window_type,pltspr)
    %This function calculates the BiSpectrum for a time series 'sample'
    %using the indirect method given the number of maximum lags 'maxshifts'
    %that are going to be used the parameters K and M for the spliting of 
    %the original sample to subsamples and the window type that is going to
    %be used during the estimation. The 'window_type' variable takes string
    %values. Specifically, it can take the values 'Rectangular' or 'Parzen'
    %to degine the usage of a rectanguar or a parzen window repsectively 
    %during the estimation. If not defined automatically takes the string
    %value 'Rectangular'. The 'pltspr' variable is a logical variable which 
    %enables the display of the plots for abosloute value of the BiSpectrum
    %estimation and the distribution of the points in the complex field for
    %each subsample if true. If false or not defined (automatically takes 
    %the value false) and disables the display of the mentioned plots.
    
    if nargin == 5
        pltspr = false;
    elseif nargin == 4
        pltspr = false;
        window_type = 'Rectangular';
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

    % Step 5
    freq1 = -0.5:1/(2*maxshifts):0.5;
    freq2 = -0.5:1/(2*maxshifts):0.5;
    C3x_f1_f2 = NaN([length(freq1) length(freq2)]);
    
    switch window_type
        case 'Rectangular'
            window = RectangularWindow(shiftings);
        case 'Parzen'
            window = ParzenWindow(shiftings);
    end
    
    for f1=1:length(freq1)
        for f2=1:length(freq2)
            C3xf1f2 = 0;
            for t1=1:length(shiftings)
                for t2=1:length(shiftings)
                    C3xf1f2 = C3xf1f2+c3x_taf1_taf2(t1,t2)*window(t1,t2)*exp(-j*(2*pi*freq1(f1)*shiftings(t1)+2*pi*freq2(f2)*shiftings(t2)));
                end
            end
            C3x_f1_f2(f1,f2) = C3xf1f2;
        end
    end
    
    switch pltspr
        case true
            figure()
            surf(freq1,freq2,abs(C3x_f1_f2));
            title('AbsoluteValue(C^x_3(f_1,f_2))','Indirect Method: '+string(window_type)+' Window - K='+string(K)+', M='+string(M)+', L_3='+string(maxshifts));
            xlabel('f_1');
            ylabel('f_2');
            zlabel('BiSpectrum');
            colorbar;
            
            figure()
            plot(C3x_f1_f2,'o','color','red');
            xline(0);
            yline(0);
            title('C^x_3(f_1,f_2) - Distribution of Complex Points','Indirect Method: '+string(window_type)+' Window - K='+string(K)+', M='+string(M)+', L_3='+string(maxshifts));
            ylabel('Imaginary (Im)');
            xlabel('Real (Re)');
    end
    
    bispec = C3x_f1_f2;
end