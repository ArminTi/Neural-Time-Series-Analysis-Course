% Create a family of complex Morlet wavelets, ranging in frequencies from 2 Hz to 30 Hz in five steps.
% Convolve each wavelet with EEG data from all electrodes and from one trial.
clc
clear
load("sampleEEGdata.mat");
num_wave = 5;
f = linspace(2,30,num_wave);

srate = EEG.srate;
time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate

family_wave = zeros(num_wave, length(time));
for freq = 1:num_wave
    s     = 6/(2*pi*f(freq));
    family_wave(freq,:) = exp(2*pi*1i*f(freq).*time) .* exp(-time.^2./(2*s^2));
end

plot(time, real(family_wave(2,:)));
%% EEG data


time = EEG.times;
data = EEG.data(:,:,1);
n_conv = length(family_wave) + EEG.pnts - 1;
half_of_wavelet_size = (length(family_wave)-1)/2;


freq_data = zeros(64,n_conv);
freq_famil_wave= zeros(num_wave, n_conv);
for i = 1:num_wave
    freq_famil_wave(i,:) = fft(family_wave(i,:),n_conv);
end
for i = 1:64
    freq_data(i,:) = fft(data(i,:),n_conv);
end
%multiply = freq_famil_wave.*freq_data;

Result = zeros(64,n_conv,num_wave);

for i = 1:64
    for j = 1:num_wave
        Result(i,:,j) = ifft(freq_famil_wave(j,:).*freq_data(i,:));
        convolution_result_fft(i,:,j) = Result(i,half_of_wavelet_size+1:end-half_of_wavelet_size,j);
    end
end

% plot for comparison
figure
subplot(311)
plot(EEG.times,real(convolution_result_fft(1,:,2)))
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
title([ 'Projection onto real axis is filtered signal at ' num2str(f(2)) ' Hz.' ])

subplot(312)
plot(EEG.times,abs(convolution_result_fft(1,:,2)).^2)
xlabel('Time (ms)'), ylabel('Power (\muV^2)')
title([ 'Magnitude of projection vector squared is power at ' num2str(f(2)) ' Hz.' ])

subplot(313)
plot(EEG.times,angle(convolution_result_fft(1,:,2)))
xlabel('Time (ms)'), ylabel('Phase angle (rad.)')
title([ 'Angle of vector is phase angle time series at ' num2str(f(2)) ' Hz.' ])

%% Now lets store phase/power to the intial matrix

%  Extract power and phase from the result of complex wavelet convolution and store in a 
% time × frequency × electrodes × power/phase matrix (thus, a 640 × 5 × 64 × 2 matrix).


time_freq_elec_feat = zeros(length(time), num_wave, 64, 2);

% Loop over channels and frequencies
for chan = 1:64
    for freq = 1:num_wave
        % Extract convolution result
        comp_signal = convolution_result_fft(chan,:,freq);
        
        % Power (magnitude squared)
        time_freq_elec_feat(:, freq, chan, 1) = abs(comp_signal).^2;
        
        % Phase (angle)
        time_freq_elec_feat(:, freq, chan, 2) = angle(comp_signal);
    end
end



%% Topographical Plots

% 4. Make topographical plots of power and phase at 180 ms at all frequencies 
% (hint: you may need to use the squeeze function to remove singleton dimensions). 
% Arrange the plots in one figure with five columns for frequency and two rows for power/phase. 
% Put labels in the plot so it is clear which topographical maps correspond to which frequencies.

[~, t_idx] = min(abs(EEG.times - 180));  % find closest time index

figure;
for freq = 1:num_wave
    subplot(2, num_wave, freq)
    topoplot(squeeze(time_freq_elec_feat(t_idx, freq,:, 1)), EEG.chanlocs);
    title([num2str(f(freq)) ' Hz Power'])
    subplot(2, num_wave, freq + num_wave)
    topoplot(squeeze(time_freq_elec_feat(t_idx, freq,:, 2)), EEG.chanlocs);
    title([num2str(f(freq)) ' Hz Phase'])
end


