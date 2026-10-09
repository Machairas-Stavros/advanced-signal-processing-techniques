clf
close all
clear
clc

load('male_voice.mat')
load('female_voice.mat')

male = {};
female = {};

male_window_coeffs = [3.3;4.1;4;4.3;3.7];
female_window_coeffs = [3.4;3.55;4.65;3.89;4.84];
male_filter_coeffs = [3;3;3;3;3];
female_filter_coeffs = [3;3;3;3;3];

for i=1:length(male_voice)
    [cepstrum,fundamental_period,chosen_signal] = calculateCepstrumPreProcessing(male_voice{i}.Audio,male_voice{i}.SamplingRate,male_window_coeffs(i),"Male",male_voice{i}.Label,true);
    impulse_response_estimation = calculateInverseCepstrum(cepstrum,fundamental_period,male_filter_coeffs(i),male_voice{i}.SamplingRate,"Male",male_voice{i}.Label,true);

    male{i} = struct('Label',male_voice{i}.Label,'Signal',chosen_signal,'FundamentalPerdiod',fundamental_period/male_voice{i}.SamplingRate,'Cepstrum',cepstrum,'ImpulseResponse',impulse_response_estimation);
end

for i=1:length(female_voice)
    [cepstrum,fundamental_period,chosen_signal] = calculateCepstrumPreProcessing(female_voice{i}.Audio,female_voice{i}.SamplingRate,female_window_coeffs(i),"Female",female_voice{i}.Label,true);
    impulse_response_estimation = calculateInverseCepstrum(cepstrum,fundamental_period,female_filter_coeffs(i),female_voice{i}.SamplingRate,"Female",female_voice{i}.Label,true);

    female{i} = struct('Label',female_voice{i}.Label,'Signal',chosen_signal,'FundamentalPerdiod',fundamental_period/female_voice{i}.SamplingRate,'Cepstrum',cepstrum,'ImpulseResponse',impulse_response_estimation);
end

clear cepstrum chosen_signal fundamental_period i impulse_response_estimation window

