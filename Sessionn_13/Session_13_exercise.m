% Pick one electrode and one time segment and compute Granger prediction between that 
% electrode (the “ seed ”) and all other electrodes in that time segment. Before selecting a 
% time segment, examine the ERP from that electrode and choose a time window that, based on the ERP, 
% is likely to contain stationary data. Justify your selection of time segment and model order. 
% Show the results in a topographical map and comment on any striking or salient features you observe.
clear
clc

load sampleEEGdata.mat
%% 

seed = 'poz';
seed_dataERP = mean(EEG.data(strcmpi(seed,{EEG.chanlocs.labels}), : , :),3);
eegdata_noERP = EEG.data(strcmpi(seed,{EEG.chanlocs.labels}), : , :) - seed_dataERP;
figure
plot(EEG.times, seed_dataERP);
hold on
plot(EEG.times, eegdata_noERP(:,:,3));

% since we have an ERP in this channel; there are two options. first,
% include a time point where data seems satationary (i.e., after 600 ms).
% or subtract the ERP from the data 
%% Granger Prediction
seed = 'poz';
% Granger prediction parameters
timewin = 300; % in ms
order   =  45; % in ms
Seed_idx = find(strcmpi(seed,{EEG.chanlocs.labels}));
% temporal down-sample results (but not data!)
times2save = 700:10:1200; % in ms

% convert parameters to indices
timewin_points = round(timewin/(1000/EEG.srate));
order_points   = round(order/(1000/EEG.srate));

% convert requested times to indices
times2saveidx = dsearchn(EEG.times',times2save');

% initialize
[x2y,y2x] = deal(zeros(64,length(times2save))); % the function deal assigns inputs to all outputs
bic = zeros(64, length(times2save),15); % Bayes info criteria (hard-coded to order=15)

for ch = 1:64
    eegdata = EEG.data([Seed_idx ch],:,:);
    for timei=1:length(times2save)
         tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        % detrend and zscore all data
    for triali=1:size(tempdata,3)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
    end
%     % reshape tempdata for armorf
     tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
     % fit AR models (model estimation from bsmart toolbox)
     [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
     [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
     [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
     %     % time-domain causal estimate
     y2x(ch, timei)=log(Ex/E(1,1));
     x2y(ch, timei)=log(Ey/E(2,2));
    end
end

%% 



