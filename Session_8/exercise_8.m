load sampleEEGdata.mat


frequency2plot = 10;  % in Hz
timewin        = 400; % in ms, for stFFT
times2save     = -300:50:1000; % in ms
timewinidx = round(timewin/(1000/EEG.srate));
% create hann taper
hann_win = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end
% define frequencies
frex = linspace(0,EEG.srate/2,floor(timewinidx/2)+1);

[junk,freq2plotidx]=min(abs(frex-frequency2plot));

% initialize ITPC output matrix
itpc = zeros(length(frex),length(times2save), length(EEG.chanlocs));



% loop over time points and perform FFT
for i = 1:length(EEG.chanlocs)
    for timepointi=1:length(times2save)
    tempdat = squeeze(EEG.data(i,times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2),:));
    taperdat = tempdat.*repmat(hann_win',1,EEG.trials);
    
    fdat = fft(taperdat,[],1)/timewinidx; % 3rd input is to make sure fft is over time
    itpc(:,timepointi,i) = abs(mean(exp(1i*angle(fdat(1:floor(timewinidx/2)+1,:))),2)); % average over trials
    end  
end
%%
figure
contourf(times2save,frex,itpc(:,:,12),40,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-200 1000])
% title([ 'ITPC at sensor ' chan2use ])

figure
plot(times2save,mean(itpc(freq2plotidx-2:freq2plotidx+2,:,1),1))
title([ 'ITPC at sensor '  ', ' num2str(frequency2plot) ' Hz' ])
set(gca,'xlim',[-200 1000])