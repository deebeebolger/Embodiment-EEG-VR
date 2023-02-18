function maxonset = EEGVR_detect_audonsets(StimIn,time,fs_eeg,dirbase,subjnum,photodtrigs,stim_path)
% Date : March 2018                         Programmed by: Deirdre Bolger
% Function to detect the onset of the auditory stimulus files recorded on
% the ERGO channel. 
% Input data: StimIn - ERGO channel data (EEG.dataERGO(2,:)) (vector)
%             time - time in ms (vector) 
%             fs_eeg - sampling rate of the input data (same as eeg
%             sampling rate) (scalar)
%             dirbase - path to data for current subject (string array).
%             subjnum - title of the current subject (e.g. s07C) (string
%             array).
%             photodtrigs - onsets and offset times of photodiode signal
%             (vector).
%             stim_path - path to the auditory stimuli *.wav files (string
%             array).
% -------------------------------------------------------------------------
%% Set the paths to the stimulation files and the sound-files (*.wav files).

stimfile = strcat(subjnum,'_stimdata.txt');
stimdir = strcat(dirbase,'StimData/',stimfile);

% Open stimdata of current subject to get order of verbs according to
% Unity.
fid = [];
fid = fopen(stimdir);
condcurr = textscan(fid, '%d%s%s%s%s%s%s','CommentStyle','//');   % This should be a 1 x 6 cell array.
fclose(fid);

trialnum = length(condcurr{1,1}); % Need to extract this information from the log files output by Unity.

condnoms = condcurr{1,2};
condunique = unique(condnoms,'stable'); % Discard repetitions of conditions.

%% ENSURE THAT CERTAIN AUDITORY STIMULI ARE PROCESSED LAST.
% As the *.wav were originally recorded at a sampling rate of 44.1kHz but,
% when recorded by Biosemi were downsampled to 1024Hz, certain auditory
% stimuli have lost their identifying features and may be confused with
% other auditory stimuli when a cross-correlation is carried out.
% Thus, these problematic stimuli are processed last to reduce the search
% space when individual auditory stimuli are being recognised via
% cross-correlation.

X = find(ismember(condunique,{'secouer' 'cacher' 'lacher' 'frotter' 'pousser'}));

for i = 1:length(X)
    
    condunique{X(i)} = [];
    
end
condunique = cat(1,condunique{:},{'pousser';'lacher';'cacher';'secouer';'frotter'});

%% DETREND THE AUDITORY STIMULI DATA CHANNEL
% To make sure it has zero mean.

StimIn =detrend(StimIn,0);

%% EXTRACT GO AND NO-GO TRIALS FROM THE PHOTODIODE TRIGGERS DATA

gotrials = {photodtrigs{:,1}};

a = {'go' 'nogo'};
goindx = find(strcmp(gotrials,a{1,1}));            % find indices of "go" trials. 
nogoindx = find(strcmp(gotrials,a{1,2}));          % find indices of "no-go" trials. 
allindx = sort([goindx nogoindx]','ascend');       % concatenate both trial-types and indices in ascending order. 
trialtime = cell2mat({photodtrigs{allindx,2}});    % the photodiode onset times in seconds
trialpnts = cell2mat({photodtrigs{allindx,4}})';   % the photodiode onset in pnts

%% LOAD IN *.WAV FILES AND CREATE A CELL ARRAY WITH A CELL FOR EACH AUDITORY FILE.
% Need to resample the loaded *.wav files so that their sampling rate
% corresponds to that of the presented auditory stimulus files. This is
% required so that the cross-correlation can function. 

sndall = cell(length(condunique),1);   % Prepare cell array. 

for len = 1:length(condunique)
    
    fwav = dir(strcat(stim_path,'*.wav'));
    fnom_wav = strcat(condunique{len,1},'.wav');
    allwavs = {fwav.name};
    
    pathnom = strcat(stim_path,fnom_wav);
    [wavin,fs] = audioread(pathnom);
    
    time_stim = [0:(length(wavin)-1)].*(1/fs);
    [wavin_rs,~] = resample(wavin,time_stim,fs_eeg);
    
    sndall{len,1} = wavin_rs;
    
end

%% FOR EACH AUDITORY STIMULUS PRESENTED (ON ERGO) DETECTING THE ONSET TIME AND THE CORRECT VERB BY CROSS-CORRELATION.
% Carried out cross-correlation between the current auditory stimulus (on
% ERGO) and each loaded *.wav file. The *.wav yielding max correlation is
% used to identify the verb presented and the lag of max correlation
% corresponds to the onset time. 

triallen = 3;
pntnum = triallen*fs_eeg;  % 3 seconds in sample points
onsettimes = zeros(length(sndall),length(trialpnts));
maxcorr = zeros(length(sndall),1);
maxonset = cell(length(trialpnts),4);

for tcnt = 1:length(trialpnts)
    disp(tcnt);
    
    stimtemp = zeros(length(StimIn),1);
    stimtemp(trialpnts(tcnt):trialpnts(tcnt)+pntnum) = 1;
    stim_curr = StimIn.*stimtemp';
    
    i = zeros(length(trialpnts),1);
    
    for condcnt = 1:length(sndall)
        
        sndcurr = sndall{condcnt,1};    % current load *.wav file. 
        [wavcorr,lags] = xcorr(stim_curr,sndcurr*-1);
        [maxcorr(condcnt),i(condcnt)] = max(wavcorr);
        onsettimes(condcnt,tcnt) = time(lags(i(condcnt)));
        
    end
    
    [~,i2] = max(maxcorr);
    maxonset{tcnt,1} = onsettimes(i2,tcnt);    % First column = onset times
    
    if ~strcmp(condunique{i2},condnoms{tcnt,1}) % Second column = condition names
        
        maxonset{tcnt,2} = condnoms{tcnt,1};
        
    else
        maxonset{tcnt,2} = condunique{i2};
    end                   
    
    maxonset{tcnt,3} = tcnt;                                    % Third column = the stimulus index
    maxonset{tcnt,4} = find(ismember(time,maxonset{tcnt,1}))';  % Fourth column = the onset latencies .
    
    marks = ones(size(time))*nan;
    marks(ismember(time, onsettimes(i2,tcnt))) = 10000;
    
    
    if tcnt == 1   % Plot the auditory stimulus with onset of one stimulus marked. 
        f1 = figure;
        subplot(2,1,1);
        plot(time./1000,StimIn);
        hold on
        stem(time./1000,marks,'or')
        subplot(2,1,2)
        plot(sndcurr)
        title(maxonset{tcnt,2});
    end

    
end







end
