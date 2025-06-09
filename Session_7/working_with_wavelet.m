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

%% Question part 2 session 7

clear all
clc
load sampleEEGdata.mat

srate   = EEG.srate;
time    = -1:1/srate:1;
chan    = 47;  
trial   = randi(EEG.trials);  

% baseline time window
baseline_time = [ -400 -100 ];
[junk,baseidx(1)] = min(abs(EEG.times-baseline_time(1)));
[junk,baseidx(2)] = min(abs(EEG.times-baseline_time(2)));

frequencies = [5 25];
nCycles     = 4;
nyquist     = srate/2;
data        = double(squeeze(EEG.data(chan,:,:)));

% Initialize result matrices
pow_wavelet = zeros(length(frequencies), EEG.pnts);
sig_wavelet = zeros(length(frequencies), EEG.pnts);
pow_hilbert = zeros(length(frequencies), EEG.pnts);
sig_hilbert = zeros(length(frequencies), EEG.pnts);

for fi = 1:length(frequencies)
    f = frequencies(fi);

    % Morlet wavelet convolution
    s = nCycles / (2*pi*f);
    wavelet = exp(2*1i*pi*f.*time) .* exp(-time.^2./(2*s^2));
    n_conv = EEG.pnts*EEG.trials + length(wavelet) - 1;
    
    data_cat = reshape(data,1,[]);
    conv_res = ifft(fft(wavelet,n_conv).*fft(data_cat,n_conv));
    
    conv_res = conv_res(floor(length(wavelet)/2)+1:end-floor(length(wavelet)/2));
    conv_res = reshape(conv_res, EEG.pnts, EEG.trials);
    
    sig_wavelet(fi,:) = real(conv_res(:,trial));
    pow_wavelet(fi,:) = mean(abs(conv_res).^2,2);



    sig_wavelet(fi,:) = 10*log10(squeeze(real(sig_wavelet(fi,:))) ./ mean(sig_wavelet(fi,baseidx(1):baseidx(2)),2) );
    pow_wavelet(fi,:) = 10*log10(squeeze(pow_wavelet(fi, :)) ./ mean(pow_wavelet(fi, baseidx(1):baseidx(2)),2) );


    % Filter-Hilbert
    freqspread = 4;  
    transwid   = 0.15;

    % ffreqs  = [0 (1-transwid)*(f-freqspread) f-freqspread f+freqspread (1+transwid)*(f+freqspread) nyquist] / nyquist;
    ffreqs  = [ 0 (1-transwid)*(f-freqspread) (f-freqspread) (f+freqspread) (1+transwid)*(f+freqspread) nyquist ]/nyquist;
    ideal   = [0 0 1 1 0 0];
    % filter_order = 100;  
    filter_order = 3*round(EEG.srate/(f-freqspread))
    fkernel = firls(filter_order, ffreqs, ideal);

    
    data_filt = reshape(filtfilt(fkernel,1,data_cat), EEG.pnts,EEG.trials);  
    % analytic  = hilbert(data_filt);        

    sig_hilbert(fi,:) = real(data_filt(:,trial));
    pow_hilbert(fi,:) = mean(abs(data_filt).^2,2);

    sig_hilbert(fi,:) = 10*log10(squeeze(real(sig_hilbert(fi,:))) ./ mean(sig_hilbert(fi,baseidx(1):baseidx(2)),2) );
    pow_hilbert(fi,:) = 10*log10(squeeze(pow_hilbert(fi, :)) ./ mean(pow_hilbert(fi, baseidx(1):baseidx(2)),2) );
end

% Plot single trial
figure
for fi = 1:2
    subplot(2,2,fi)
    plot(EEG.times, sig_wavelet(fi,:), 'k'); hold on
    plot(EEG.times, sig_hilbert(fi,:), 'r')
    set(gca,'xlim',[-300 800])
    title(['Filtered signal at ' num2str(frequencies(fi)) ' Hz'])
    legend('Wavelet','Hilbert')
    xlabel('Time (ms)'), ylabel('Amplitude')
    grid on

    subplot(2,2,fi+2)
    plot(EEG.times, pow_wavelet(fi,:), 'k'); hold on
    plot(EEG.times, pow_hilbert(fi,:), 'r')
    set(gca,'xlim',[-300 800])
    title(['Power at ' num2str(frequencies(fi)) ' Hz'])
    legend('Wavelet','Hilbert')
    xlabel('Time (ms)'), ylabel('Power (μV^2)')
    grid on
end
