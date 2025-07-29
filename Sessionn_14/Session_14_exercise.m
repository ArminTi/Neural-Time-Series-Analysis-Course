% 1. Pick two electrodes and two frequencies (one frequency per electrode) and compute
% mutual information over time and trials between power from the first electrode and phase
% from the second electrode. Justify your choice of bin size. Next, recompute mutual information using phase from the first electrode and power from the second electrode. Make sure you
% use the same bin size you used in the previous analysis, so the results are directly comparable.
% Plot the time courses of the mutual information from these two analyses.

clear
clc

load sampleEEGdata.mat

% Channel of interest
channel_1 = 'Fz';
channel_2 = 'C5';
electrodesidx(1) = find(strcmpi(channel_1,{EEG.chanlocs.labels}));
electrodesidx(2) = find(strcmpi(channel_2,{EEG.chanlocs.labels}));

% Frequency of interest
Freq_1 = 5;
Freq_2 = 13;

% specify convolution and wavelet info
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

% FFT of data
fft_EEG1 = fft(reshape(EEG.data(electrodesidx(1),:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_EEG2 = fft(reshape(EEG.data(electrodesidx(2),:,:),1,EEG.pnts*EEG.trials),n_convolution);

% FFT of wavelet
fft_wavelet1 = fft(exp(2*1i*pi*Freq_1.*time) .* exp(-time.^2./(2*(4/(2*pi*Freq_1))^2)),n_convolution);
fft_wavelet2 = fft(exp(2*1i*pi*Freq_2.*time) .* exp(-time.^2./(2*(4/(2*pi*Freq_2))^2)),n_convolution);

% Specify time
timewindow = 300; % in ms
times2save = -300:20:800;
timewindowidx = round(timewindow/(1000/EEG.srate)/2);
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end
%% Fz 5 hz power and C5 13 hz phase

convres     = ifft(fft_wavelet1.*fft_EEG1,n_convolution);
fz_5   = reshape(convres(half_of_wavelet_size+1:end-half_of_wavelet_size),EEG.pnts,EEG.trials);
convres     = ifft(fft_wavelet2.*fft_EEG2,n_convolution);
C5_13   = reshape(convres(half_of_wavelet_size+1:end-half_of_wavelet_size),EEG.pnts,EEG.trials);


% initialize outputs
mi   = zeros(2, length(times2save)); % Fz 5 hz power and C5 13 hz phase


for timei = 1:length(times2save)
        datax = fz_5(times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx,:);
        datay = C5_13(times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx,:);

        power = log10(abs(datax).^2);
        phase = angle(datay);
        
        % compute MI
        mi(1, timei) = mutualinformationx(power,phase,50); %fixed bin
        mi(2, timei) = mutualinformationx(power,phase); %variable bin
end




%% Plotting

subplot(121)
plot(times2save, mi(1,:));
xlabel('Time (ms)'), ylabel('MI (bits)')
title('Constant bin length')
set(gca,'xlim',[-300 800],'ylim',[min(mi(:))-.01 max(mi(:))+.01])

subplot(122)
plot(times2save, mi(2,:));
xlabel('Time (ms)'), ylabel('MI (bits)')
title('Variable bin length')
set(gca,'xlim',[-300 800],'ylim',[min(mi(:))-.01 max(mi(:))+.01])

