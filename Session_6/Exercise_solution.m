
% For question part 6
clear
clc
load sampleEEGdata.mat

% Wave and data let info
srate = 500;
time = -1:1/srate:1;
lowest_frequency  =   2;  % in Hz
highest_frequency = 30;  % in Hz
frequencies=linspace(lowest_frequency,highest_frequency,5);


% concate data
data = EEG.data(1,:,:);
data2filter_cat = squeeze(double(reshape(data,1,EEG.pnts*EEG.trials)));
nConv = length(data2filter_cat) + length(time) - 1;
fft_data = fft(data2filter_cat, nConv);


tf_triali_final = zeros(length(frequencies), EEG.pnts, EEG.trials);
% tf_avg = zeros(length(frequencies), num_cov);
for fi = 1:length(frequencies)
    wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*(6/(2*pi*frequencies(fi)))^2));
    wave_fft = fft(wavelet, nConv);
    tf_triali = ifft(fft_data.*wave_fft);
    tf_triali = real(tf_triali);  
    % Trim edges
    half_wave = floor(length(time)/2);
    tf_triali = tf_triali(half_wave+1:end-half_wave);
    
    % Reshape back to [time x trials]
    tf_triali = reshape(tf_triali, EEG.pnts, EEG.trials);
    tf_triali_final(fi,:,:) = tf_triali;
end
%% Plotting


% Average the result of convolution over all trials and plot an ERP corresponding to each wavelet frequency. 
% Each frequency should be in its own subplot.

Average_tf_final = mean(tf_triali_final(:,:,:),3);

for i = 1:length(frequencies)
    subplot(5,1,i);
    plot(EEG.times,Average_tf_final(i, :));
end
%% Comparing with the ERP
subplot(6,1,1);
plot(EEG.times,mean(data,3));

subplot(6,1,2);  
plot(EEG.times,Average_tf_final(1, :));

subplot(6,1,3);  
plot(EEG.times,Average_tf_final(2, :));

subplot(6,1,4);  
plot(EEG.times,Average_tf_final(3, :));

subplot(6,1,5);  
plot(EEG.times,Average_tf_final(4, :));

subplot(6,1,6);  
plot(EEG.times,Average_tf_final(5, :));











