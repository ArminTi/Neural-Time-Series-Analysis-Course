clear
clc

% loading
load sampleEEGdata.mat
sampling_rate = EEG.srate;
channel_1 = 'Fz';
channel_2 = 'O1';
Freqs = [5, 10, 15];

% Time vector
timevec = EEG.times / 1000; % convert to seconds
n_times = length(timevec);
n_trials = EEG.trials;

% Extract data & reshape
chan1_idx = find(strcmpi(channel_1, {EEG.chanlocs.labels}));
chan2_idx = find(strcmpi(channel_2, {EEG.chanlocs.labels}));
data1 = reshape(EEG.data(chan1_idx,:,:), EEG.pnts, n_trials);
data2 = reshape(EEG.data(chan2_idx,:,:), EEG.pnts, n_trials);
times2save = -300:20:800; % in ms
times2saveidx = dsearchn(EEG.times',times2save');

% Parameters
wavelet_cycles = 4.5;
time = -1:1/sampling_rate:1;
half_wavelet = (length(time)-1)/2;
n_wavelet = length(time);
n_data = EEG.pnts * EEG.trials;
n_conv = n_wavelet + n_data - 1;

fft_data1 = fft(reshape(data1,1,n_data), n_conv);
fft_data2 = fft(reshape(data2,1,n_data), n_conv);

% Define window types
win_types = {'3cycles', '150ms', '900ms'};


tf_corrdata_all = zeros(3,length(times2save),length(win_types));

for fi = 1:length(Freqs)
    freq = Freqs(fi);

    % Create wavelet
    s = wavelet_cycles / (2*pi*freq);
    wavelet = exp(2*1i*pi*freq.*time) .* exp(-time.^2/(2*s^2));
    fft_wavelet = fft(wavelet, n_conv);

    % Convolve & extract power 
    conv_result1 = ifft(fft_wavelet .* fft_data1, n_conv) * sqrt(wavelet_cycles /(2*pi*freq));
    conv_result1 = conv_result1(half_wavelet+1:end-half_wavelet);
    conv_result1 = abs(reshape(conv_result1, EEG.pnts, n_trials)).^2;

    conv_result2 = ifft(fft_wavelet .* fft_data2, n_conv) * sqrt(wavelet_cycles /(2*pi*freq));
    conv_result2 = conv_result2(half_wavelet+1:end-half_wavelet);
    conv_result2 = abs(reshape(conv_result2, EEG.pnts, n_trials)).^2;

    for w = 1:length(win_types)
        % Compute window length in points
        switch win_types{w}
            case '3cycles'
                win_ms = 3*(1000/freq); %win_ms = 1000/freq; % I dont know why, but it seems cohen did not multiply by three in its solution figures. yet, it seems a bit wrong!
            case '150ms'
                win_ms = 150;
            case '900ms'
                win_ms = 900;
        end

        timewin_points = round(win_ms/(1000/sampling_rate));

        for timei = 1:length(times2save)
            idx1 = times2saveidx(timei)-floor(timewin_points/2);
            idx2 = times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2);

            tempdata_1 = mean(conv_result1(idx1:idx2,:),1);
            tempdata_2 = mean(conv_result2(idx1:idx2,:),1);

            r = corr(tempdata_1', tempdata_2','type', 'Spearman');
            tf_corrdata_all(fi,timei,w) = r;
        end
    end
end

%% Plotting
figure

colors = {'b','g','r'};
win_labels = {'150 ms','900 ms','3f'};

for fi = 1:length(Freqs)
    subplot(3,1,fi)
    hold on
    for w = 1:length(win_types)
        % For 3f, compute actual window length at this frequency
        if strcmp(win_types{w},'3cycles')
            win_len = round(3*(1000/Freqs(fi))); % in ms
            label_str = sprintf('3f (%d ms)', win_len);
        else
            label_str = sprintf('%s', win_types{w});
        end
        
        plot(times2save, tf_corrdata_all(fi,:,w), ...
            'Color',colors{w},'LineWidth',1.5,'DisplayName',label_str)
    end
    set(gca,'xlim',[-300,800])
    xlabel('Time (ms)')
    ylabel('Spearman \rho')
    ylim([-0.3 0.3])
    title(sprintf('Correlation between %s and %s at %d Hz',channel_1,channel_2,Freqs(fi)))
    legend
    grid on
end

sgtitle('Time-resolved power correlation at different frequencies and window sizes')

