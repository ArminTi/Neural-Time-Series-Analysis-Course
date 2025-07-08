% 1. Select one seed electrode and one frequency band and compute phase-based connectivity between that seed electrode 
% and every other electrode. Use two methods for phase-based connectivity that were presented in this chapter, 
% one that is volume conduction independent (e.g., PLI) and one that could produce spurious connectivity 
% due to volume conduction (e.g., ISPC). Do not apply a baseline subtraction. Make topographical plots 
% of seeded connectivity in a time window of your choice (e.g., 300–350 ms). 
% What are the similarities and differences between results from the two methods, and 
% what might be the reasons for the similarities and differences?  
clear
clc
load sampleEEGdata.mat

% names of the channels you want to synchronize
channel1 = 'Cz';
chanidx_seed = find(strcmpi(channel1,{EEG.chanlocs.labels}));

% create complex Morlet wavelet
center_freq = 8; % in Hz
time        = -1:1/EEG.srate:1; % time for wavelet
wavelet     = exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(4/(2*pi*center_freq))^2))/center_freq;
half_of_wavelet_size = (length(time)-1)/2;

% time of interest
times2save = -400:20:800;
times2saveidx = dsearchn(EEG.times',times2save');

% FFT parameters
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

% FFT of wavelet
fft_wavelet = fft(wavelet,n_convolution);

% initialize output time-frequency data
phase_data = zeros(64,EEG.pnts);
real_data  = zeros(64,EEG.pnts);


% fft for seed
fft_seed = fft(reshape(EEG.data(chanidx_seed,:,:),1,n_data),n_convolution);
convolution_result_fft_seed = ifft(fft_wavelet.*fft_seed,n_convolution) * sqrt(4/(2*pi*center_freq));
convolution_result_fft_seed = convolution_result_fft_seed(half_of_wavelet_size+1:end-half_of_wavelet_size);
phase_data_seed = angle(reshape(convolution_result_fft_seed,EEG.pnts,EEG.trials));

% run convolution and extract filtered signal (real part) and phase

for chani=1:64
    fft_data = fft(reshape(EEG.data(chani,:,:),1,n_data),n_convolution);
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(4/(2*pi*center_freq));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    phase_data = angle(reshape(convolution_result_fft,EEG.pnts,EEG.trials));

    % phase angle differences
    phase_diffs = phase_data-phase_data_seed;
    
    % compute ICPS over time
    ps(chani,:) = abs(mean(exp(1i*phase_diffs(times2saveidx,:)),2));

end
%% Topoplot


timewin = [300 350];
timeidx = times2save >= timewin(1) & times2save <= timewin(2);

% Average ICPS over the chosen time window
icps_avg = mean(ps(:,timeidx),2);

figure;
topoplot(icps_avg, EEG.chanlocs, ...
    'maplimits', [0 1], ...
    'electrodes', 'on');
colorbar;
title(sprintf('Seeded Connectivity with %s (%d–%d ms)', channel1, timewin(1), timewin(2)));



