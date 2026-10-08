clf
close all
clear 
clc

rng('default');

N = 2048;
p = 0;
d = 0;
q = 5;
q_sub = q-2;
q_sup = q+3;
ma_coeffs = [1,0.93,0.85,0.72,0.59,-0.10];
mean_noise = 1;
K = 32;
M = 64;
L3 = 20;
shiftings = -L3:L3;
SNR = 30:-5:-5;
no_iterations = 50;

signal_x = zeros(N,no_iterations);
signal_v = NaN(N,no_iterations);
skewness_v = NaN(1,no_iterations);
c3x_t1_t2 = NaN(length(shiftings),length(shiftings),no_iterations);
impulse_response_x_q = NaN(q+1,no_iterations);
signal_x_estim_q = NaN(N,no_iterations);
nrmse_x_q = NaN(1,no_iterations);
impulse_response_x_qsub = NaN(q_sub+1,no_iterations);
signal_x_estim_qsub = NaN(N,no_iterations);
nrmse_x_qsub = NaN(1,no_iterations);
impulse_response_x_qsup = NaN(q_sup+1,no_iterations);
signal_x_estim_qsup = NaN(N,no_iterations);
nrmse_x_qsup = NaN(1,no_iterations);
signal_y = NaN(size(signal_x,1),length(SNR),no_iterations);
c3y_t1_t2 = NaN(size(c3x_t1_t2,1),size(c3x_t1_t2,2),length(SNR),no_iterations);
impulse_response_y =NaN(q+1,length(SNR),no_iterations);
nrmse_y = NaN(length(SNR),no_iterations);
signal_y_estim = NaN(N,length(SNR),no_iterations);

for iter=1:no_iterations
    fprintf('%.f\n',iter)
    
    signal_v(:,iter) = exprnd(mean_noise,[N,1]);
    
    for i=1:N
        for j=0:q
            if i<=j
                break
            end
            signal_x(i,iter) = signal_x(i,iter)+ma_coeffs(j+1)*signal_v(i-j,iter);
        end
    end

    skewness_v(iter) = sum((signal_v(:,iter)-mean(signal_v(:,iter))).^3)/((length(signal_v(:,iter))-1)*std(signal_v(:,iter))^3);

    [c3x_t1_t2(:,:,iter),impulse_response_x_q(:,iter),signal_x_estim_q(:,iter),nrmse_x_q(iter)] = cumulantsImpulseNRMSE(signal_x(:,iter),signal_v(:,iter),q,K,M,L3);
    [~,impulse_response_x_qsub(:,iter),signal_x_estim_qsub(:,iter),nrmse_x_qsub(iter)] = cumulantsImpulseNRMSE(signal_x(:,iter),signal_v(:,iter),q_sub,K,M,L3);
    [~,impulse_response_x_qsup(:,iter),signal_x_estim_qsup(:,iter),nrmse_x_qsup(iter)] = cumulantsImpulseNRMSE(signal_x(:,iter),signal_v(:,iter),q_sup,K,M,L3);    

    for i=1:length(SNR)
        signal_y(:,i,iter) = awgn(signal_x(:,iter),SNR(i),'measured');
        [c3y_t1_t2(:,:,i,iter),impulse_response_y(:,i,iter),signal_y_estim(:,i,iter),nrmse_y(i,iter)] = cumulantsImpulseNRMSE(signal_y(:,i,iter),signal_v(:,iter),q,K,M,L3);
    end

end

mean_signal_x = mean(signal_x,2);
mean_signal_v = mean(signal_v,2);
mean_skewness_v = mean(skewness_v);
mean_c3x_t1_t2 = mean(c3x_t1_t2,3);
mean_nrmse_x_q = mean(nrmse_x_q);
mean_impulse_response_x_q = mean(impulse_response_x_q,2);
mean_signal_x_estim_q = mean(signal_x_estim_q,2);
mean_nrmse_x_qsub = mean(nrmse_x_qsub);
mean_impulse_response_x_qsub = mean(impulse_response_x_qsub,2);
mean_signal_x_estim_qsub = mean(signal_x_estim_qsub,2);
mean_nrmse_x_qsup = mean(nrmse_x_qsup);
mean_impulse_response_x_qsup = mean(impulse_response_x_qsup,2);
mean_signal_x_estim_qsup = mean(signal_x_estim_qsup,2);
mean_signal_y = mean(signal_y,3);
mean_c3y_t1_t2 = mean(c3y_t1_t2,4);
mean_nrmse_y = mean(nrmse_y,2);
mean_impulse_response_y = mean(impulse_response_y,3);
mean_signal_y_estim = mean(signal_y_estim,3);

figure()
tiledlayout(2,1)
nexttile()
plot(1:N,mean_signal_x);
title('Mean Real Discrete Signal x[k] - MA Process - Number Of Iterations: '+string(no_iterations));
nexttile()
plot(1:N,mean_signal_v);
title('Mean White Non-Gaussian Noise v[k] - Noise - Number Of Iterations: '+string(no_iterations))

figure()
surf(shiftings,shiftings,mean_c3x_t1_t2);
title('Mean 3^{rd} Order Cumulants - Number Of Iterations: '+string(no_iterations),'Indirect Method: K='+string(K)+', M='+string(M)+', L_3='+string(L3));
xlabel('\tau_1');
ylabel('\tau_2');
zlabel('c^x_3(\tau_1,\tau_2)');
colorbar;

figure()
tiledlayout(3,1)
nexttile()
stem(1:q+1,mean_impulse_response_x_q)
title('Mean Impulse Response - q = '+string(q));
subtitle('Real Value Of q');
nexttile()
stem(1:q_sub+1,mean_impulse_response_x_qsub)
title('Mean Impulse Response - q = '+string(q_sub))
subtitle('Sub-Estimation Of Order q - Deviation From Real Value: '+string(q-q_sub));
nexttile()
stem(1:q_sup+1,mean_impulse_response_x_qsup)
title('Mean Impulse Response - q = '+string(q_sup))
subtitle('Sup-Estimation Of Order q - Deviation From Real Value: '+string(q-q_sup)');
sgtitle('\bf{Mean Impulse Response In Function Of Order q - Without Additive Noise - Number Of Iterations: }'+string(no_iterations));

figure()
tiledlayout(3,1)
nexttile()
plot(1:N,mean_signal_x,'color','blue');
hold on
plot(1:N,mean_signal_x_estim_q,'color','red');
text('String','MeanNRMSE='+string(mean_nrmse_x_q),'FontSize',9,'Position',[2050,mean(mean_signal_x)]);
hold off
title('Mean x[k] vs x_{est}[k]');
legend('x[k]','x_{est}[k]');
nexttile()
plot(1:N,mean_signal_x,'color','blue');
hold on;
plot(1:N,mean_signal_x_estim_qsub,'color','red');
text('String','MeanNRMSE='+string(mean_nrmse_x_qsub),'FontSize',9,'Position',[2050,mean(mean_signal_x)]);
hold off;
title('Mean x[k] Vs Mean x_{est,sub}[k]','q_{sub}='+string(q_sub));
legend('x[k]','x_{est,sub}[k]');
nexttile()
plot(1:N,mean_signal_x,'color','blue');
hold on;
plot(1:N,mean_signal_x_estim_qsup,'color','red');
text('String','MeanNRMSE='+string(mean_nrmse_x_qsup),'FontSize',9,'Position',[2050,mean(mean_signal_x)]);
hold off;
title('Mean x[k] Vs Mean x_{est,sup}[k]','q_{sup}='+string(q_sup));
legend('x[k]','x_{est,sup}[k]');
sgtitle('\bf{Mean Estimations Of x[k] For Real Value & Sub/Sup-Estimation Of q}')

figure()
for i=1:length(SNR)
    nexttile()
    plot(1:N,mean_signal_y(:,i));
    title('SNR='+string(SNR(i))+'dB');
end
sgtitle('\bf{Mean Output Signal y_i[k]=x[k]+n_i[k] - Additive White Gaussian Noise - Number Of Iterations: }'+string(no_iterations));

figure()
for i=1:length(SNR)
    nexttile()
    surf(shiftings,shiftings,mean_c3y_t1_t2(:,:,i));
    title('Indirect Method: K='+string(K)+', M='+string(M)+', L_3='+string(L3));
    subtitle('SNR='+string(SNR(i))+'dB');
    xlabel('\tau_1');
    ylabel('\tau_2');
    zlabel('c^y_3(\tau_1,\tau_2)');
    colorbar;
end
sgtitle('\bf{Mean 3^{rd} Order Cumulants C^{y_i}_3 - Additive White Gaussian Noise - Number Of Iterations: }'+string(no_iterations));

figure()
for i=1:length(SNR)
    nexttile()
    plot(1:N,mean_signal_y(:,i),'color','blue');
    hold on;
    plot(1:N,mean_signal_y_estim(:,i),'color','red');
    text('String','MeanNRMSE='+string(mean_nrmse_y(i)),'FontSize',9,'Position',[2050,mean(mean_signal_y(:,i))]);
    hold off;
    title('SNR='+string(SNR(i))+'dB');
    legend('x[k]','x_{est}[k]');
end
sgtitle('\bf{Mean x[k] Vs x_{estim,i}[k] - Number Of Iterations: }'+string(no_iterations));

figure()
for i=1:length(SNR)
    nexttile()
    stem(1:length(mean_impulse_response_y(:,i)),mean_impulse_response_y(:,i));
    title('SNR='+string(SNR(i))+'dB');
end
sgtitle('\bf{Mean Impulse Response h_i[k] - Additive White Gaussian Noise - Number Of Iterations: }'+string(no_iterations));

figure()
plot(SNR,mean_nrmse_y,'Marker','o','color','red');
set(gca,'XDir','reverse');
title('Mean Normalized Root Mean Square Error In Function Of SNR - Number Of Iterations: '+string(no_iterations))
xlabel('SNR (dB)');
ylabel('NRMSE');

clc
fprintf('MeanSkewness: %f\n',mean_skewness_v);
fprintf('MeanNRMSE: %f - q: %.f\n',mean_nrmse_x_q,q);
fprintf('MeanNRMSE: %f - q_sub: %.f\n',mean_nrmse_x_qsub,q_sub);
fprintf('MeanNRMSE: %f - q_sup: %.f\n',mean_nrmse_x_qsup,q_sup);
fprintf('Mean MeanNRMSE for all SNR values: %f - q: %.f\n',mean(mean_nrmse_y),q);