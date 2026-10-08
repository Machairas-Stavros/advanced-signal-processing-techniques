clf;
close all;
clear;
clc;

no_iterations = 50;

powerspec = [];
bispec_indirect = [];
bispec_indiparz = [];
bispec_direct = [];
x_k = [];

for i=1:no_iterations
    fprintf('%.f\n',i);
    %1] Construct the X[k]
    lamda = NaN([1 6]);
    phase = NaN([1 6]);

    lamda(1) = 0.12;
    lamda(2) = 0.30;
    lamda(3) = lamda(1)+lamda(2);
    lamda(4) = 0.19;
    lamda(5) = 0.17;
    lamda(6) = lamda(4)+lamda(5);

    phase(1) = unifrnd(0,2*pi);
    phase(2) = unifrnd(0,2*pi);
    phase(3) = phase(1)+phase(2);
    phase(4) = unifrnd(0,2*pi);
    phase(5) = unifrnd(0,2*pi);
    phase(6) = phase(4)+phase(5);

    omega = 2*pi*lamda;

    N = 8192;
    k = 0:N-1; 
    
    endcol = size(x_k,2);
    for a=1:N
        x_k(a,endcol+1) = sum(cos(omega*k(a)+phase));
    end

    x_k = x_k-mean(x_k); 
    %The mean values is extracted in order to produce a stationary signal. The 
    %estimated cumulants are going to be expressed with a simplier equation -
    %relation.
    
    %2] Estimate the Power Spectrum C^x_2(f). Use L2=128 max shifting for the
    %autocorrelation
    L2 = 128;
    powerspec(:,end+1) = PowerSpectrum(x_k(:,end),L2);

    %3] Estimate the BiSpectrum (only in the primary area) C^3_x(f1,f2) using:
    K = 32;
    M = 256;
    L3 = 64;

    %a] The Indrirect Method with K=23,M=256. Use L3=64 max shiftings for the
    %third order cumulants 

    %a1] Use Rectangular Window
    bispec_indirect(:,:,end+1) = BiSpectrumIndirect(x_k(:,end),K,M,L3,'Rectangular');

    %a2] Use Parzen Window
    bispec_indiparz(:,:,end+1) = BiSpectrumIndirect(x_k(:,end),K,M,L3,'Parzen');

    %b] The Direct Method with K=32,M=256,J=0.
    J = 0;
    fs = 2^14;
    bispec_direct(:,:,end+1) = BiSpectrumDirect(x_k(:,end),K,M,fs,J);
end

figure()
plot(0:1:N-1,mean(x_k,2),'Marker','o');
title('Mean(x[k])', 'Number of Iterations: '+string(no_iterations));
ylabel('Discrete Process');
xlabel('k=0,...,8191');

MeanSpecs(L2,L3,K,M,J,fs,no_iterations,powerspec,bispec_indirect,bispec_indiparz,bispec_direct,true);