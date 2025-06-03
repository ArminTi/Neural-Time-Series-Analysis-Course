%% Question 1
% 1. Pick three electrodes. Compute time-frequency plots of ITPC and decibel-corrected power
% for these electrodes, using either complex Morlet wavelet convolution or the filter-Hilbert
% method. Plot the results side by side for each electrode (power and ITPC in subplots; one
% figure for each electrode). Are the patterns of results from ITPC and power generally similar
% or generally different? Do the results look more similar at some electrodes and less similar at
% other electrodes?

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


% convwavelet to fft of data in each channel
t_p = zeros(length(elec_idx), length(freq), EEG.pnts);
ITPC = zeros(length(elec_idx), length(freq), EEG.pnts);
for ch = 1:length(elec_idx)
    for f = 1:length(freq)
        wave_let = exp(2*1i*pi*freq(f).*time) .* exp(-time.^2./(2*(4/(2*pi*freq(f)))^2))/freq(f);
        fft_wave = fft(wave_let,n_conv);
        convolution_result_fft = ifft(fft_wave.*data_fft(ch,:),n_conv);
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
        t_p(ch,f,:) = mean(abs(convolution_result_fft).^2,2);

        % base norm
        t_p(ch,f,:) = 10*log10(squeeze(t_p(ch,f,:) ./ mean(t_p(ch,f,baseidx(1):baseidx(2)),3)));

        % Averaging phase (in literature you might see calling this
        % coherence)
        ITPC(ch,f,:) = abs(mean(exp(1i*angle(convolution_result_fft)),2));
    end
end

% Plotting
analysis_labels = {'Cz', 'Fp1', 'O1'};
clims = repmat([-3 3], length(elec_idx), 1); 
figure;
for ch = 1:length(elec_idx)
    % Plot raw power
    subplot(length(elec_idx), 2, (ch-1)*2 + 1);
    contourf(EEG.times, freq, squeeze(t_p(ch,:,:)), 40, 'linecolor', 'none');
    set(gca, 'clim', clims(ch,:), 'xlim', [-400 1000], 'xtick', -200:200:800);
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    title(['Raw Power - ' analysis_labels{ch}]);
    colorbar;

    % Plot ITPC
    subplot(length(elec_idx), 2, (ch-1)*2 + 2);
    contourf(EEG.times, freq, squeeze(ITPC(ch,:,:)), 40, 'linecolor', 'none');
    set(gca, 'xlim', [-400 1000], 'xtick', -200:200:800);
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    title(['ITPC - ' analysis_labels{ch}]);
    colorbar;
end

sgtitle('Time-Frequency Analysis: Raw Power & ITPC');



%% Question 2
% 2. For each of these three electrodes, compute wITPCz using reaction time as the trial-varying
% modulator. Perform this analysis for all time-frequency points to generate time-frequency
% maps of the relationship between phase and reaction time. Do the time-frequency maps
% of wITPCz look different from the time-frequency maps of ITPC? Do you see any striking
% patterns in the ITPCz results, and do the results differ across the different electrodes (don’t
% worry about statistics, base your judgment on qualitative patterns)? How would you interpret
% the results if they were statistically significant?