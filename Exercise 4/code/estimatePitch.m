function pitch_estim = estimatePitch(impresp,voicesig)
    
    impresp = vertcat(impresp,zeros(length(impresp),1));
    
    voicesig_fft = fft(voicesig);
    voicesig_fft = fftshift(voicesig_fft);

    impresp_fft = fft(impresp);
    impresp_fft = fftshift(impresp_fft);

    faketrans_function = voicesig_fft./impresp_fft;

    pitch_estim = ifft(faketrans_function);
    pitch_estim = ifftshift(pitch_estim);

end