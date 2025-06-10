% Perform a time-frequency decomposition of the data from two electrodes both before and after computing 
% the surface Laplacian (that is, compute the surface Laplacian on the raw data before applying a 
% time-frequency decomposition). Compute both power (decibels normalized using a baseline period of your choice) 
% Plot the results using the same color scaling for before and after the surface Laplacian.
% Are there any salient differences in the time-frequency power or ITPC results before versus after application of the surface Laplacian, and do the differences depend on the frequency?
% How would you interpret similarities and differences at different frequency bands?  
clear
clc
% Load data from three electrod
load sampleEEGdata.mat
electrode_of_inter = {'Cz', 'Fp1', 'O1'};
all_labels = {EEG.chanlocs.labels};
elec_idx = find(ismember(lower(all_labels), lower(electrode_of_inter)));
data = EEG.data(elec_idx, :, :);

% wavelet info
low_f = 2;
high_f = 30;
num_f = 20;
freq = logspace(log10(low_f),log10(high_f),num_f);
srate = EEG.srate;
time = -1:1/srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% baseline
baseline = [-400 -100];
[junk,baseidx(1)] = min(abs(EEG.times-baseline(1))); % convert baseline from ms to indices
[junk,baseidx(2)] = min(abs(EEG.times-baseline(2)));

% FFT of data
n_wave = length(time);
n_data = EEG.pnts*EEG.trials;
n_conv =  n_wave + n_data -1;
for ch = 1:length(elec_idx)
    % reshape single channel across trials into 1D
    temp_data = reshape(squeeze(data(ch, :, :)), 1, []);
    data_fft(ch, :) = fft(temp_data, n_conv);
end

% Laplace of data

X = [EEG.chanlocs.X];
Y = [EEG.chanlocs.Y];
Z = [EEG.chanlocs.Z];
data_lap = laplacian_perrinX(EEG.data,X,Y,Z);

for ch = 1:length(elec_idx)
    % reshape single channel across trials into 1D
    temp_data = reshape(squeeze(data_lap(ch, :, :)), 1, []);
    data_lap_fft(ch, :) = fft(temp_data, n_conv);
end



%% convwavelet to fft of data in each channel
t_p = zeros(length(elec_idx), length(freq), EEG.pnts);
ITPC = zeros(length(elec_idx), length(freq), EEG.pnts);

t_p_lap = zeros(length(elec_idx), length(freq), EEG.pnts);
ITPC_lap = zeros(length(elec_idx), length(freq), EEG.pnts);
for ch = 1:length(elec_idx)
    for f = 1:length(freq)
        % Wavelet
        wave_let = exp(2*1i*pi*freq(f).*time) .* exp(-time.^2./(2*(4/(2*pi*freq(f)))^2))/freq(f);
        fft_wave = fft(wave_let,n_conv);

        % conv with raw data
        convolution_result_fft = ifft(fft_wave.*data_fft(ch,:),n_conv);
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
        t_p(ch,f,:) = mean(abs(convolution_result_fft).^2,2);
        t_p(ch,f,:) = 10*log10(squeeze(t_p(ch,f,:) ./ mean(t_p(ch,f,baseidx(1):baseidx(2)),3)));
        ITPC(ch,f,:) = abs(mean(exp(1i*angle(convolution_result_fft)),2));

        % conv with lap data
        convolution_result_fft = ifft(fft_wave.*data_lap_fft(ch,:),n_conv);
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
        t_p_lap(ch,f,:) = mean(abs(convolution_result_fft).^2,2);
        t_p_lap(ch,f,:) = 10*log10(squeeze(t_p_lap(ch,f,:) ./ mean(t_p_lap(ch,f,baseidx(1):baseidx(2)),3)));
        ITPC_lap(ch,f,:) = abs(mean(exp(1i*angle(convolution_result_fft)),2));


    end
end
%% Plotting

figure(1);

clims = repmat([-3 3], length(elec_idx), 1); 
subplot(321)
contourf(EEG.times, freq, squeeze(t_p(1,:,:)), 40, 'linecolor', 'none');
set(gca, 'clim', clims(ch,:), 'xlim', [-400 1000], 'xtick', -200:200:800);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Raw Power - Before Laplace Cz');
colorbar;

subplot(322)
contourf(EEG.times, freq, squeeze(t_p_lap(1,:,:)), 40, 'linecolor', 'none');
set(gca, 'clim', clims(ch,:), 'xlim', [-400 1000], 'xtick', -200:200:800);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Raw Power - After Laplace Cz');
colorbar;

subplot(323)
contourf(EEG.times, freq, squeeze(t_p(2,:,:)), 40, 'linecolor', 'none');
set(gca, 'clim', clims(ch,:), 'xlim', [-400 1000], 'xtick', -200:200:800);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Raw Power - Before Laplace Fp1');
colorbar;

subplot(324)
contourf(EEG.times, freq, squeeze(t_p_lap(2,:,:)), 40, 'linecolor', 'none');
set(gca, 'clim', clims(ch,:), 'xlim', [-400 1000], 'xtick', -200:200:800);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Raw Power - After Laplace Fp1');
colorbar;

subplot(325)
contourf(EEG.times, freq, squeeze(t_p(3,:,:)), 40, 'linecolor', 'none');
set(gca, 'clim', clims(ch,:), 'xlim', [-400 1000], 'xtick', -200:200:800);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Raw Power - Before Laplace O1');
colorbar;

subplot(326)
contourf(EEG.times, freq, squeeze(t_p_lap(3,:,:)), 40, 'linecolor', 'none');
set(gca, 'clim', clims(ch,:), 'xlim', [-400 1000], 'xtick', -200:200:800);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Raw Power - After Laplace O1');
colorbar;






