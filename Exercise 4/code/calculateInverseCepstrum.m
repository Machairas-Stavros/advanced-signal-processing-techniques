function impresp_estim = calculateInverseCepstrum(cepstrum,fundamental_period,noperiods_winlen,samplingRate,gender,label,pltspr)
    
    if nargin==5
        pltspr = false;
    end

    window_length = round(noperiods_winlen*fundamental_period);
    window_center = floor(length(cepstrum)/2)+1;
     if mod(window_length,2)==0
        window_start = window_center-floor(window_length/2)+1;
        window_end = window_center+floor(window_length/2);
    else
        window_start = window_center-floor(window_length/2);
        window_end = window_center+floor(window_length/2);
     end
    window = zeros(length(cepstrum),1);
    window(window_start:window_end) = rectwin(window_length);

    liftered_cepstrum = cepstrum.*window;
    liftered_cepstrum_cut = liftered_cepstrum(liftered_cepstrum~=0);

    signal_to_fft = liftered_cepstrum(floor(length(liftered_cepstrum)/2)+1:end);
    nozerospadded = 2^nextpow2(length(signal_to_fft))-length(signal_to_fft);
    signal_to_fft = vertcat(signal_to_fft,zeros(nozerospadded,1));
    liftered_cepstrum_fft = fft(signal_to_fft);
    liftered_cepstrum_fft = fftshift(liftered_cepstrum_fft);
    liftered_cepstrum_exp_fft = power(10,liftered_cepstrum_fft);

    impresp_estim = ifft(liftered_cepstrum_exp_fft);
    impresp_estim = ifftshift(impresp_estim);

    freq_axis_liftcep = (-length(liftered_cepstrum_fft)/2:length(liftered_cepstrum_fft)/2-1)*samplingRate/length(liftered_cepstrum_fft);
    time_axis_liftcep = (-length(liftered_cepstrum)/2:length(liftered_cepstrum)/2-1)*1/samplingRate;
    time_axis_liftcep_cut = (-length(liftered_cepstrum_cut)/2:length(liftered_cepstrum_cut)/2-1)*1/samplingRate;
    time_axis_impresp_estim = (-length(impresp_estim)/2:length(impresp_estim)/2-1)*1/samplingRate;

    switch pltspr
        case true
            figure()
            nexttile()
            plot(time_axis_liftcep_cut,window(window_start:window_end))
            ylabel('Amplitude')
            xlabel('Quefrency (sec)')
            title('Rectangular Window')
            subtitle('Low Pass Filter - LTI')
            
            nexttile()
            plot(time_axis_liftcep,cepstrum)
            ylabel('Amplitude')
            xlabel('Quefrency (sec)')
            title('Complex Cepstrum - IFFT(Log_{10}(|FFT(PPOVS)|))')
        
            nexttile()
            plot(time_axis_liftcep,liftered_cepstrum)
            ylabel('Amplitude')
            xlabel('Quefrency (sec)')
            title('Liftered Complex Cepstrum - LTI(Cepstrum)')
        
            nexttile()
            plot(time_axis_liftcep_cut,liftered_cepstrum_cut)
            ylabel('Amplitude')
            xlabel('Quefrency (sec)')
            title('Liftered Complex Cepstrum Cut - LTI(|Cepstrum|)')
            subtitle('Omitting Zero Values - High Frequencies')
        
            sgtitle(gender+' - Letter "'+label+'"');
        
            figure()
            nexttile()
            hold on
            yline(0)
            xline(0)
            plot(liftered_cepstrum_fft,'o')
            hold off
            ylabel('Imaginary (Im)')
            xlabel('Real (Re)')
            title('Complex Points Distribution - FFT(LTI(Cepstrum))')
        
            nexttile()
            plot(freq_axis_liftcep,abs(liftered_cepstrum_fft))
            ylabel('Amplitude (dB)')
            xlabel('Frequency (Hz)')
            title('Absolute Value FFT(LTI(Cepstrum))')
        
            nexttile()
            plot(freq_axis_liftcep,phaseReduction(phase(liftered_cepstrum_fft)))
            yticks([-pi,-pi/2,0,pi/2,pi])
            yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
            ylabel('Phase (rad)')
            xlabel('Frequency (Hz)')
            title('Phase FFT(LTI(Cepstrum))')
        
            nexttile()
            hold on
            yline(0)
            xline(0)
            plot(liftered_cepstrum_exp_fft,'o')
            hold off
            ylabel('Imaginary (Im)')
            xlabel('Real (Re)')
            title('Complex Points Distribution - 10^{FFT(LTI(Cepstrum))}')
        
            nexttile()
            plot(freq_axis_liftcep,abs(liftered_cepstrum_exp_fft))
            ylabel('Amplitude')
            xlabel('Frequency (Hz)')
            title('Absolute Value 10^{FFT(LTI(Cepstrum))}')
        
            nexttile()
            plot(freq_axis_liftcep,phaseReduction(phase(liftered_cepstrum_exp_fft)))
            yticks([-pi,-pi/2,0,pi/2,pi])
            yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
            ylabel('Phase (rad)')
            xlabel('Frequency (Hz)')
            title('Phase 10^{FFT(LTI(Cepstrum))}')
        
            sgtitle(gender+' - Letter "'+label+'"');
        
            figure()
            nexttile()
            plot(time_axis_impresp_estim,impresp_estim)
            ylabel('Amplitude')
            xlabel('Time (sec)')
            title('Complex Value Impulse Response Estimation - IFFT(10^{FFT(LTI(Cepstrum))})')
        
            nexttile()
            plot(time_axis_impresp_estim,abs(impresp_estim))
            ylabel('Amplitude')
            xlabel('Time (sec)')
            title('Absolute Value Impulse Response Estimation - IFFT(10^{FFT(LTI(Cepstrum))})')
        
            nexttile()
            plot(time_axis_impresp_estim,phaseReduction(phase(impresp_estim)))
            yticks([-pi,-pi/2,0,pi/2,pi])
            yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
            ylabel('Phase (rad)')
            xlabel('Time (sec)')
            title('Phase Impulse Response Estimation - IFFT(10^{FFT(LTI(Cepstrum))})')
        
            sgtitle(gender+' - Letter "'+label+'"');
    end
end