clf
close all
clear 
clc

rng('default');

N = 2048;
p = 0;
d = 0;
q = 5;
ma_coeffs = [1,0.93,0.85,0.72,0.59,-0.10];

mean_noise = 1;
signal_v = exprnd(mean_noise,[N,1]);

signal_x = zeros(N,1);

for i=1:N
    for j=0:q
        if i<=j
            break
        end
        signal_x(i) = signal_x(i)+ma_coeffs(j+1)*signal_v(i-j);
    end
end

figure()
tiledlayout(2,1)
nexttile()
plot(1:N,signal_x);
title('Real Discrete Signal x[k] - MA Process');
nexttile()
plot(1:N,signal_v);
title('White Non-Gaussian Noise v[k] - Noise')

skewness_v = sum((signal_v-mean(signal_v)).^3)/((length(signal_v)-1)*std(signal_v)^3);
fprintf('Skewness: %f\n',skewness_v);

K = 32;
M = 64;
L3 = 20;

shiftings = -L3:L3;

q_sub = q-2;
q_sup = q+3;

[c3x_t1_t2,impulse_response_x_q,signal_x_estim_q,nrmse_x_q] = cumulantsImpulseNRMSE(signal_x,signal_v,q,K,M,L3);
[~,impulse_response_x_qsub,signal_x_estim_qsub,nrmse_x_qsub] = cumulantsImpulseNRMSE(signal_x,signal_v,q_sub,K,M,L3);
[~,impulse_response_x_qsup,signal_x_estim_qsup,nrmse_x_qsup] = cumulantsImpulseNRMSE(signal_x,signal_v,q_sup,K,M,L3);

fprintf('NRMSE: %f - q: %.f\n',nrmse_x_q,q);
fprintf('NRMSE: %f - q_sub: %.f\n',nrmse_x_qsub,q_sub);
fprintf('NRMSE: %f - q_sup: %.f\n',nrmse_x_qsup,q_sup);

figure()
surf(shiftings,shiftings,c3x_t1_t2);
title('3^{rd} Order Cumulants','Indirect Method: K='+string(K)+', M='+string(M)+', L_3='+string(L3));
xlabel('\tau_1');
ylabel('\tau_2');
zlabel('c^x_3(\tau_1,\tau_2)');
colorbar;

figure()
tiledlayout(3,1)
nexttile()
stem(1:q+1,impulse_response_x_q)
title('Impulse Response - q = '+string(q));
subtitle('Real Value Of q');
nexttile()
stem(1:q_sub+1,impulse_response_x_qsub)
title('Impulse Response - q = '+string(q_sub))
subtitle('Sub-Estimation Of Order q - Deviation From Real Value: '+string(q-q_sub));
nexttile()
stem(1:q_sup+1,impulse_response_x_qsup)
title('Impulse Response - q = '+string(q_sup))
subtitle('Sup-Estimation Of Order q - Deviation From Real Value: '+string(q-q_sup)');
sgtitle('\bf{Impulse Response In Function Of Order q - Without Additive Noise}');

figure()
tiledlayout(3,1)
nexttile()
plot(1:N,signal_x,'color','blue');
hold on;
plot(1:N,signal_x_estim_q,'color','red');
text('String','NRMSE='+string(nrmse_x_q),'FontSize',9,'Position',[2050,mean(signal_x)]);
hold off;
title('x[k] Vs x_{est}[k]','q='+string(q));
legend('x[k]','x_{est}[k]');
nexttile()
plot(1:N,signal_x,'color','blue');
hold on;
plot(1:N,signal_x_estim_qsub,'color','red');
text('String','NRMSE='+string(nrmse_x_qsub),'FontSize',9,'Position',[2050,mean(signal_x)]);
hold off;
title('x[k] Vs x_{est,sub}[k]','q_{sub}='+string(q_sub));
legend('x[k]','x_{est,sub}[k]');
nexttile()
plot(1:N,signal_x,'color','blue');
hold on;
plot(1:N,signal_x_estim_qsup,'color','red');
text('String','NRMSE='+string(nrmse_x_qsup),'FontSize',9,'Position',[2050,mean(signal_x)]);
hold off;
title('x[k] Vs x_{est,sup}[k]','q_{sup}='+string(q_sup));
legend('x[k]','x_{est,sup}[k]');
sgtitle('\bf{Estimations Of x[k] For Real Value & Sub/Sup-Estimation Of q}')

SNR = 30:-5:-5;

signal_y = NaN(length(signal_x),length(SNR));
c3y_t1_t2 = NaN(2*L3+1,2*L3+1,length(SNR));
impulse_response_y =NaN(q+1,length(SNR));
signal_y_estim = NaN(N,length(SNR));
nrmse_y = NaN(1,length(SNR));

for i=1:length(SNR)
    signal_y(:,i) = awgn(signal_x,SNR(i),'measured');
    [c3y_t1_t2(:,:,i),impulse_response_y(:,i),signal_y_estim(:,i),nrmse_y(i)] = cumulantsImpulseNRMSE(signal_y(:,i),signal_v,q,K,M,L3);
end

fprintf('Mean NRMSE for all SNR values: %f - q: %.f\n',mean(nrmse_y),q);

figure()
for i=1:length(SNR)
    nexttile()
    plot(1:N,signal_y(:,i));
    title('SNR='+string(SNR(i))+'dB');
end
sgtitle('\bf{Output Signal y_i[k]=x[k]+n_i[k] - Additive White Gaussian Noise}');

figure()
for i=1:length(SNR)
    nexttile()
    plot(1:N,signal_y(:,i),'color','blue');
    hold on;
    plot(1:N,signal_y_estim(:,i),'color','red');
    hold off;
    text('String','NRMSE='+string(nrmse_y(i)),'FontSize',9,'Position',[2050,mean(signal_y(:,i))]);
    title('SNR='+string(SNR(i))+'dB');
    legend('x[k]','x_{est}[k]');
end
sgtitle('\bf{x[k] Vs x_{estim,i}[k]}');

figure()
for i=1:length(SNR)
    nexttile()
    surf(shiftings,shiftings,c3y_t1_t2(:,:,i));
    title('Indirect Method: K='+string(K)+', M='+string(M)+', L_3='+string(L3));
    subtitle('SNR='+string(SNR(i))+'dB');
    xlabel('\tau_1');
    ylabel('\tau_2');
    zlabel('c^y_3(\tau_1,\tau_2)');
    colorbar;
end
sgtitle('\bf{3^{rd} Order Cumulants C^{y_i}_3 - Additive White Gaussian Noise}');

figure()
for i=1:length(SNR)
    nexttile()
    stem(1:length(impulse_response_y(:,i)),impulse_response_y(:,i));
    title('SNR='+string(SNR(i))+'dB');
end
sgtitle('\bf{Impulse Response h_i[k] - Additive White Gaussian Noise}');

figure()
plot(SNR,nrmse_y,'Marker','o','color','red');
set(gca,'XDir','reverse');
title('Normalized Root Mean Square Error In Function Of SNR')
xlabel('SNR (dB)');
ylabel('NRMSE');
