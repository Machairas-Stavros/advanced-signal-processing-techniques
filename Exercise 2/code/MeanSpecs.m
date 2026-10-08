function [mpowerspec,mbispec_indirect,mbispec_indiparz,mbispec_direct] = MeanSpecs(L2,L3,K,M,J,fs,noi,powerspec,bispec_indirect,bispec_indiparz,bispec_direct,pltspr)
    %This function takes as an input the two and three dimesnional matrices
    %derived from the stage of the 50 iterations for loop and calculates,
    %plots and returns the meann values of each matrix through the second
    %and the third dimension repsectively. The resulted new matrices
    %represent the mean value of the Power Spectrum and BiSpectrum
    %repsectively taking into account the results of each iteration. The
    %'pltspr' variable is used to enable or not the creation of the plots
    %if true or false respectively.
    
    if nargin == 11
        pltspr = false;
    end
    
    freq_ps = -0.5:1/(2*L2):0.5;
    freq_bsid = -0.5:1/(2*L3):0.5;
    freq_bsd = -0.5:1/(M-1):0.5;

    
    mpowerspec = mean(powerspec,2);
    mbispec_indirect = mean(bispec_indirect,3);
    mbispec_indiparz = mean(bispec_indiparz,3);
    mbispec_direct = mean(bispec_direct,3);

    switch pltspr
        case true
            figure()
            plot(freq_ps,abs(mpowerspec));
            title('Mean(AbsoluteValue(C^x_2(f_1,f_2)))','Number Of Iterations: '+string(noi)+' - L_2='+string(L2));
            xlabel('frequency');
            ylabel('Power Spectrum');
            
            figure()
            surf(freq_bsid,freq_bsid,abs(mbispec_indirect));
            title('Mean(AbsoluteValue(C^x_3(f_1,f_2)))','Number Of Iterations: '+string(noi)+' - Indirect Method: Rectangular Window, K='+string(K)+', M='+string(M)+', L_3='+string(L3));
            xlabel('f_1');
            ylabel('f_2');
            zlabel('BiSpectrum');
            colorbar;
            
            figure()
            surf(freq_bsid,freq_bsid,abs(mbispec_indiparz));
            title('Mean(AbsoluteValue(C^x_3(f_1,f_2)))','Number Of Iterations: '+string(noi)+' - Indirect Method: Parzen Window, K='+string(K)+', M='+string(M)+', L_3='+string(L3));
            xlabel('f_1');
            ylabel('f_2');
            zlabel('BiSpectrum');
            colorbar;
            
            figure()
            surf(freq_bsd,freq_bsd,abs(mbispec_direct));
            title('Mean(AbsoluteValue(C^x_3(f_1,f_2)))','Number Of Iterations: '+string(noi)+' - Direct Method, K='+string(K)+', M='+string(M)+', J='+string(J)+', f_s='+string(fs));
            xlabel('f_1');
            ylabel('f_2');
            zlabel('BiSpectrum');
            colorbar;
    end
end