function [cepstrum,fundamental_period,chosen_signal] = calculateCepstrumPreProcessing(ogsig,samplingRate,noperiods_winlen,gender,label,pltspr)
    
    if nargin==5
        pltspr=false;
    end

    prepro_ogsig = ogsig-mean(ogsig);
    prepro_ogsig = prepro_ogsig/abs(max(prepro_ogsig));

    autocorr_prepro_ogsig = autocorr(prepro_ogsig,'NumLags',length(prepro_ogsig)-1);
    [autocorr_pks,autocorr_pks_locs] = findpeaks(autocorr_prepro_ogsig);
    [~,ind_maxpeak] = max(autocorr_pks);
    fundamental_period = autocorr_pks_locs(ind_maxpeak);
    fundamental_period = fundamental_period-1;
    fundamental_period_sec = fundamental_period/samplingRate;
    fundamental_frequency_hz = 1/fundamental_period_sec;

    [~,window_center] = max(prepro_ogsig);
    window_length = round(noperiods_winlen*fundamental_period);
    if mod(window_length,2)==0
        window_start = window_center-window_length/2+1;
        window_end = window_center+window_length/2;
    else
        window_start = window_center-floor(window_length/2);
        window_end = window_center+floor(window_length/2);
    end
    window = hamming(window_length);

    chosen_signal = prepro_ogsig(window_start:window_end);
    windowed_signal = chosen_signal.*window;

    nozerospadded = 2^nextpow2(length(windowed_signal))-length(windowed_signal);
    ppovs = vertcat(windowed_signal,zeros(nozerospadded,1));
    
    ppovs_fft = fft(ppovs);
    ppovs_fft = fftshift(ppovs_fft);    
    ppovs_log_fft = log10(abs(ppovs_fft));
    cepstrum = ifft(ppovs_log_fft);
    cepstrum = ifftshift(cepstrum);

    freq_axis = (-length(ppovs_fft)/2:length(ppovs_fft)/2-1)*samplingRate/length(ppovs_fft);
    time_axis_prepro_ogsig = (0:length(ogsig)-1)*(1/samplingRate);
    time_axis_cep = (-length(ppovs_fft)/2:length(ppovs_fft)/2-1)*1/samplingRate;
    time_axis_window = (window_start:window_end)*1/samplingRate;
    time_axis_ppovs = (window_start:window_end+nozerospadded)*1/samplingRate;

    switch pltspr
        case true
            figure()
            nexttile()
            plot(time_axis_prepro_ogsig,prepro_ogsig)
            ylabel('Amplitude')
            xlabel('Time (sec)')
            title('Original Voice Signal')
            subtitle('Normalized & Mean Extracted')
        
            nexttile()
            plot(time_axis_window,chosen_signal)
            ylabel('Amplitude')
            xlabel('Time (sec)')
            title('Choosen Part Of The Signal')
            subtitle('Fundamental Period: '+string(fundamental_period_sec)+'sec')
        
            nexttile()
            plot(time_axis_window,hamming(window_length))
            ylabel('Amplitude')
            xlabel('Time (sec)')
            title('Hamming Window')
        
            nexttile()
            plot(time_axis_window,windowed_signal)
            ylabel('Amplitude')
            xlabel('Time (sec)')
            title('Windowed Signal')
        
            nexttile()
            plot(time_axis_ppovs,ppovs)
            ylabel('Amplitude')
            xlabel('Time (sec)')
            title('Pre-Processed Original Voice Signal (PPOVS)')
            subtitle('Zero-Padded Windowed Signal')

            sgtitle(gender+' - Letter "'+label+'"')
        
            figure()
            nexttile()
            hold on
            xline(0)
            yline(0)
            plot(ppovs_fft,'o')
            hold off
            ylabel('Imaginary (Im)')
            xlabel('Real (Re)')
            title('Complex Points Distribution - FFT(PPOVS)')
        
            nexttile()
            plot(freq_axis,abs(ppovs_fft))
            ylabel('Amplitude')
            xlabel('Frequency (Hz)')
            title('Absolute Value - |FFT(PPOVS)|')
        
            nexttile()
            plot(freq_axis,phaseReduction(phase(ppovs_fft)))
            yticks([-pi,-pi/2,0,pi/2,pi])
            yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
            ylabel('Phase (rad)')
            xlabel('Frequency (Hz)')
            title('Phase - FFT(PPOVS)')
        
            nexttile()
            plot(freq_axis,ppovs_log_fft)
            ylabel('Amplitude (dB)')
            xlabel('Frequency (Hz)')
            title('Common Logarithm - Log_{10}(|FFT(PPOVS)|)')
        
            sgtitle(gender+' - Letter "'+label+'"')
    
            figure()
            nexttile()
            plot(time_axis_cep,cepstrum)
            ylabel('Amplitude')
            xlabel('Quefrency (sec)')
            title('Complex Cepstrum - IFFT(Log_{10}(|FFT(PPOVS)|))')
        
            nexttile()
            plot(time_axis_cep,abs(cepstrum))
            ylabel('Amplitude')
            xlabel('Quefrency (sec)')
            title('Absolute Value Of Cepstrum - |IFFT(Log_{10}(|FFT(PPOVS)|))|')
        
            nexttile()
            plot(time_axis_cep,phaseReduction(phase(cepstrum)))
            yticks([-pi,-pi/2,0,pi/2,pi])
            yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
            ylabel('Phase (rad)')
            xlabel('Quefrency (sec)')
            title('Phase Of Cepstrum - IFFT(Log_{10}(|FFT(PPOVS)|))')
        
            sgtitle(gender+' - Letter "'+label+'"')
    end
end