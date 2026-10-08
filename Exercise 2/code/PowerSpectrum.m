function powerspec = PowerSpectrum(sample,maxshifts,pltspr)
    %This function calculates the Power Spectrum for a time series 'sample'
    %given the number of maximum lags 'maxshifts' that are going to be
    %used. the 'pltspr' variable is a logical variable which enables the 
    %display of the plots for the autocorellation, the distribution of the 
    %points in the complex field, the absolute value of the Power Spectrum
    %and the phase of the Power Spectrum if true. If it has the value false
    %or if it is not defined (automatically takes the value false) and
    %disables the display of the mentioned plots.
    %pltspr takes the values true or the value false
    if nargin == 2
        pltspr = false;
    end
    
    shiftings = -maxshifts:1:maxshifts;
    
    c2x_taf = NaN([1 length(shiftings)]);
    for i=1:length(shiftings)
        c2x_taf(i) = (1/length(sample))*sum(sample.*arrayShift(sample,shiftings(i)));    
    end

    f = -0.5:1/(2*maxshifts):0.5;
    C2x_f = NaN([1 length(f)]);
    for i=1:length(f)
        C2x_f(i) = sum(c2x_taf.*exp(-j*2*pi*f(i)*shiftings));
    end
    
    switch pltspr
        case true
            figure()
            stem(shiftings,c2x_taf);
            title('c^x_2(\tau)','AutoCorrelation - L_2='+string(maxshifts));
            ylabel('AutoCorrelation');
            xlabel('Shiftings');
            
            figure()
            plot(f,abs(C2x_f));
            title('AbsoluteValue(C^x_2(f))','L_2='+string(maxshifts));
            ylabel('Power Spectrum');
            xlabel('Frequency');
            
            figure()
            plot(f,phase(C2x_f));
            title('Phase(C^x_2(f))','L_2='+string(maxshifts));
            ylabel('Power Spectrum');
            xlabel('Frequency');

            figure()
            plot(C2x_f,'o','color','red');
            xline(0);
            yline(0);
            title('C^x_2(f) - Distribution of Complex Points','L_2='+string(maxshifts));
            ylabel('Imaginary (Im)');
            xlabel('Real (Re)');
    end
    
    powerspec = C2x_f;
end