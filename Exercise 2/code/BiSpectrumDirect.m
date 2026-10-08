function bispec = BiSpectrumDirect(sample,K,M,fs,J,pltspr)
    %This function calculates the BiSpectrum for a time series 'sample'
    %using the direct method given the parameters K and M for the spliting
    %of the original sample to subsamples and the J and fs parameters for 
    %the smoothing that is going to be used during the estimation. J is 
    %equal to zero, in this case, so no smoothing is implemented. The 
    %'pltspr' variable is a logical variable which enables the display of 
    %the plots for the abosloute value of the BiSpectrum estimation and the
    %distribution of the points in the complex field for each subsample if 
    %true. If false or not defined (automatically takes the value false) 
    %and disables the display of the mentioned plots.
    
    if nargin == 5
        pltspr = false;
    elseif nargin == 4
        pltspr = false;
        J = 0;
    elseif nargin == 3
        pltspr = false;
        J = 0;
        fs = 1;
    end
    
    D = 2*J+1;
    N0 = round(M/D);
    delta0 = fs/N0;
    
    m = 0:1:M-1;
    
    % Step 1
    subsamples = NaN([M K]);
    for i=0:K-1
        subsamples(:,i+1) = sample(i*M+1:(i+1)*M);
    end

    % Step 2
    subsamples = subsamples-mean(subsamples);

    % Step 3
    dft_subsamples = NaN([M K]);
    for i=1:K
        for lamda=0:M-1
            dft_subsamples(lamda+1,i) = sum(subsamples(:,i).'.*exp(-j*2*pi*lamda*m/M));
        end
    end

    % Step 4: Smoothing is not implemented in the code

    % Step 5
    b3x_l1_l2 = zeros([M M K]);
    for k=1:K
        for l1=1:M
            for l2=1:M
                b3x_l1_l2(l1,l2,k) = (1/delta0^2)*dft_subsamples(l1,k)*dft_subsamples(l2,k)*conj(dft_subsamples(mod(l1+l2,M)+1,k));
            end
        end
    end

    % Step 6
    C3x_w1_w2 = fftshift(sum(b3x_l1_l2,3));

    switch pltspr
        case true
            ff = -0.5:1/(M-1):0.5;

            figure()
            surf(ff,ff,abs(C3x_w1_w2));
            title('AbsoluteValue(C^x_3(f_1,f_2))','Direct Method - K='+string(K)+', M='+string(M)+', J='+string(J)+', f_s='+string(fs));
            xlabel('f_1');
            ylabel('f_2');
            zlabel('BiSpectrum');
            colorbar;
            
            figure()
            plot(C3x_w1_w2,'o','color','red');
            xline(0);
            yline(0);
            title('C^x_3(w_1,w_2) - Distribution of Complex Points','Direct Method - K='+string(K)+', M='+string(M)+', J='+string(J)+', f_s='+string(fs));
            ylabel('Imaginary (Im)');
            xlabel('Real (Re)');
    end
    
    bispec = C3x_w1_w2;
end