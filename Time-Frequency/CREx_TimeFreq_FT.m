% Script to carry out time frequency analysis on input EEG data...
% 02-2018.   Programmed by: D. Bolger
% This script uses fieldtrip functions to carry out the time-frequency decomposition.
% It consists of two parts: a first which prepares the segmented data for time-frequency decomposition.
%                           a second part in which time-frequency decomposition is carried out.
%**************************************************************

clear all
close all
AllChans=load(''); % Load the channel info file.

TFcfg=[];

TFcfg.currdir= fullfile(filesep,);  % Path to data 
TFcfg.Condnom1 ={};                 % Enter the conditions
TFcfg.Condnom2 ={};
TFcfg.Condnom3 ={};
TFcfg.Subnum = ;    % Enter the number of subjects                                                                                                                                      % Number of subjects.
TFcfg.Channum = 64; % Enter the number of electrodes
TFcfg.extend_mirror='no';
TFcfg.maxfreq = 80;   % Define maximum frequency of interest
trialextend='no'; 

%Create the full directory path
fullpaths=cell(length(TFcfg.Condnom1),1);
condnames=cell(length(TFcfg.Condnom1),1);

for cnt=1:length(TFcfg.Condnom1)
    if ~isempty(TFcfg.Condnom3)
        fullpaths{cnt,1}=fullfile(TFcfg.currdir,TFcfg.Condnom1{1,cnt},TFcfg.Condnom2{1,cnt},TFcfg.Condnom3{1,cnt},filesep);
        condnames{cnt,1}=strcat(TFcfg.Condnom1{1,cnt},'-',TFcfg.Condnom2{1,cnt},'-',TFcfg.Condnom3{1,cnt});
    elseif ~isempty(TFcfg.Condnom2)
        fullpaths{cnt,1}=fullfile(TFcfg.currdir,TFcfg.Condnom1{1,cnt},TFcfg.Condnom2{1,cnt},filesep);
        condnames{cnt,1}=strcat(TFcfg.Condnom1{1,cnt},'-',TFcfg.Condnom2{1,cnt});
    elseif isempty(TFcfg.Condnom2)
        fullpaths{cnt,1}=fullfile(TFcfg.currdir,TFcfg.Condnom1{1,cnt},filesep);
        condnames{cnt,1}=strcat(TFcfg.Condnom1{1,cnt});
    end
end

%% Load the data using EEGLAB interface

files_curr = cell(length(TFcfg.Condnom1),1);
fileIndx   = cell(length(TFcfg.Condnom1),1);
AllData    = cell(TFcfg.Subnum,length(TFcfg.Condnom1));
AllData_FT = cell(TFcfg.Subnum,length(TFcfg.Condnom1));
AllData_FText = cell(TFcfg.Subnum,length(TFcfg.Condnom1));
tdata_extend  = cell(TFcfg.Subnum,length(TFcfg.Condnom1));
tdata_in      = cell(TFcfg.Subnum,length(TFcfg.Condnom1));
T             =cell(TFcfg.Subnum,length(TFcfg.Condnom1));
padlimit   = 5000;  %define padding limits in ms.
targ_event = {};

for condcnt=1:length(TFcfg.Condnom1)
    
    [ALLEEG, ~, ~, ALLCOM] = eeglab;  %open an EEGLAB session
    
    display(fullpaths{condcnt,1});                          %display the contents of the current folder
    files_curr{condcnt,1} = dir(fullpaths{condcnt,1});
    currfiles = files_curr{condcnt,1};
    fileIndx{condcnt,1} = find(~[currfiles.isdir]);
    FInd = fileIndx{condcnt,1};
    filenum = dir(strcat(fullpaths{condcnt,1},'*.set'));
    filenom = {filenum.name};                                 %titles of the all the .set files to be loaded
    
    EEG = pop_loadset('filename',filenom,'filepath',fullpaths{condcnt,1});
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
    EEG=eeg_checkset(EEG); eeglab redraw;
    
    for counter=1:length(ALLEEG)
        
        AllData{counter,condcnt}    = ALLEEG(counter);
        AllData_FT{counter,condcnt} = CREx_EEGFT_configure(ALLEEG(counter),fullpaths{condcnt,1},filenom{1,counter},'no',targ_event); % Call of function to configure EEGLAB data for use in FieldTrip
        
        if strcmp(trialextend,'yes')
            tdata = cell2mat(AllData_FT{counter,condcnt}.trial);
            tdata_in{counter,condcnt} = reshape(tdata,[length(AllData_FT{1,1}.label) size(tdata,2)/length(AllData_FT{counter,condcnt}.trial) length(AllData_FT{counter,condcnt}.trial)]);
            
            [tdata_extend{counter,condcnt}, tdata_lims,T{counter,condcnt}]= CREx_trialextend_mirror(tdata_in{counter,condcnt},'both',AllData_FT{counter,condcnt}.fsample,AllData_FT{counter,condcnt}.time{1,1});
            
            X = squeeze(mat2cell(tdata_extend{counter,condcnt},size(tdata_extend{counter,condcnt},1),size(tdata_extend{counter,condcnt},2),...
                ones(1,size(tdata_extend{counter,condcnt},3))))';
            T_ext = repmat(T{counter,condcnt},[1,1,size(tdata_extend{counter,condcnt},3)]);
            XT = squeeze(mat2cell(T_ext,1,size(T{counter,condcnt},2),ones(1,size(tdata_extend{counter,condcnt},3))))';
            sampleinfo=zeros(size(tdata_extend{counter,condcnt},3),2);
            col1=1;
            for xcnt=1:size(sampleinfo,1)
                sampleinfo(xcnt,1)= col1;
                sampleinfo(xcnt,2)= (size(tdata_extend{counter,condcnt},2)*xcnt)+xcnt-1;
                col1=sampleinfo(xcnt,2)+2;
            end
            FTdata_new = [];
            FTdata_new = AllData_FT{counter,condcnt};
            FTdata_new.sampleinfo = sampleinfo;
            FTdata_new.trial = X;
            FTdata_new.time = XT;
            FTdata_new.hdr.nSamples = size(tdata_extend{counter,condcnt},2);
            FTdata_new.hdr.nSamplesPre = find([T{counter,condcnt}==0]);
            AllData_FText{counter,condcnt} = FTdata_new;
        else
            AllData_FText{counter,condcnt}=AllData_FT{counter,condcnt};
        end
    end
end
clear condcnt counter

%% Calculate the Grand Average over all trials for each individual subject
% The mean over all trials can be used to calculate the evoked
% Time-Frequency.

tlcfg=[];
tlcfg.channel = {AllChans.chaninfo(1:72).labels};
tlcfg.trials='all';
tlcfg.keeptrials='no';

% Initialise the data arrays.
AllData_mean = cell(TFcfg.Subnum,length(TFcfg.Condnom1));
AllData_induced = AllData_FT;
colnum = size(AllData_FT{1,1}.trial{1,1},2);

for condcnt=1: length(TFcfg.Condnom1)
    for counter=1:TFcfg.Subnum
        [AllData_mean{counter,condcnt}] = ft_timelockanalysis(tlcfg,AllData_FT{counter,condcnt});   %Call of function to calculate the GA for each subject
        AllData_mean{counter,condcnt}.fsample = AllData_FT{counter,condcnt}.fsample;
        
        colnum = size(AllData_FT{counter,condcnt}.trial{1,1},2);
        trialnum = size(AllData_FT{counter,condcnt}.trial,2);
        Dtrials = reshape(cell2mat(AllData_FT{counter,condcnt}.trial),[],colnum,trialnum);
        Devoked = repmat(AllData_mean{counter,condcnt}.avg,1,1,trialnum);
        Induced_pow = Dtrials - Devoked;
        v= squeeze(num2cell(Induced_pow,[1 2]))';
        AllData_induced{counter,condcnt}.trial = v;
    end
end


%% Determine the frequencies that can be resolved given the current trial length and sampling frequency 

trialtime=AllData_FT{1,1}.time{1,1};
SR=AllData_FT{1,1}.fsample; 


freqpos = 0:(1/trialtime(end)):(SR/2); 
disp(['Current sampling frequency is ',num2str(SR),'Hz']) 
disp(['Current trial length is ', num2str(trialtime(end)),'seconds'])
disp(['The lowest frequency that you can resolve is ',num2str(freqpos(2)*2),'Hz']); 

pad1 = ceil(numel(trialtime)/SR);
txtra=pad1*(1/SR);
padder=(trialtime(end)+abs(trialtime(1))+txtra)*4;   %Calculate zero padding... 

nfft=2^nextpow2((numel(trialtime)+pad1)*2);
[freqs,~]=freqgrid_calc(SR,nfft,[freqpos(2) TFcfg.maxfreq],10);       %Call of function to define the frequency grid

%%  Choose TF-decomposition method and define necessary parameters
disp('-------Time to carry out the time-frequency decomposition-------------'); 
methodchoice={'MultiTaper' 'Wavelet Convolution' 'CREx Complex Morlet Wavelet'};
disp('Choose one of the following time-frequency decomposition methods: ');
disp('MultiTaper (FieldTrip) , Morlet Wavelet Convolution (FieldTrip), Morlet Wavelet (CREx):')
choice=input('Choose "1" for MultiTaper, "2" for Wavelet Convolution or "3" for Crex Complex MorletWavelet :');

disp(['*****Carrying out Time-Frequency decomposition using ',methodchoice{1,choice},' method*****']);

if strcmp(methodchoice{1,choice},'MultiTaper')
    TFout_mt = MultiTaper_FT(AllData_FT,freqs,padder);                      % Use of FieldTrip function to carry out TF decomposition using multi-tapers
elseif strcmp(methodchoice{1,choice},'Wavelet Convolution')                  % Use of FielTrip function to carry out TF decomposition using Morlet wavelet convolution
     TFout_mixed = MWaveletConv_FT(AllData_FT,freqs,padder,'all');              % Calculate the evoked+induced time-frequency
     TFout_evoked = MWaveletConv_FT(AllData_mean,freqs,padder,'mean');            % Calculate the evoked time-frequency       
     TFout_induced = MWaveletConv_FT(AllData_induced,freqs,padder,'all');
elseif strcmp(methodchoice{1,choice},'CREx Complex Morlet Wavelet')
    %CREx_waveletcomplex()                                                   %Use of home-made TF calculation using complex Morlet wavelets.
end

%% CALCULATE INDUCED BY SUBTRACTING THE MEAN TIME-FREQUENCY DATA (EVOKED) FROM SINGLE-TRIAL TF DATA
TFout_inducedv2 = TFout_mixed;
TFout_avg_inducedv2 = cell(size(TFout_mixed));
TFGA_avg_inducedv2 = cell(1,size(TFout_mixed,2));
twindow = [-0.2 1.0];

fprintf('**Calculating the induced activity by subtracting evoked tf data from single-trial tf data**');

fdindcfg = [];
fdindcfg.keeptrials = 'no';    % find average across trials
fdindcfg.channels = 'all';
fdindcfg.frequency = 'all';
fdindcfg.latency = twindow;

GAcfg = [];
GAcfg.keepindividual = 'yes';
GAcfg.foilim = 'all';
GAcfg.toilim = twindow;
GAcfg.channel = cellstr(TFout_evoked{1,1}.label(1:72));
GAcfg.parameter = 'powspctrm';

for condcnt = 1:size(TFout_mixed,2)
    fprintf('***condition number %d\n',condcnt);
    for sujcnt = 1:size(TFout_mixed,1)
        fprintf('***subject %d\n',sujcnt)
        tnum_curr = size(TFout_mixed{sujcnt,condcnt}.powspctrm,1);
        tosub = repmat(TFout_evoked{sujcnt,condcnt}.powspctrm,1,1,1,tnum_curr);
        tosub = permute(tosub,[4 1 2 3]);
        TFout_inducedv2{sujcnt,condcnt}.powspctrm = TFout_mixed{sujcnt,condcnt}.powspctrm - tosub;
        TFout_avg_inducedv2{sujcnt,condcnt} = ft_freqdescriptives(fdindcfg,TFout_inducedv2{sujcnt,condcnt});
    end
    TFGA_avg_inducedv2{1,condcnt} = ft_freqgrandaverage(GAcfg,TFout_avg_inducedv2{:,condcnt}); 
    
end

%% Find the subject-level Grand Average (average over trials)
% Here the baseline is the average baseline over all trials and it is
% applied to the Grand Average time-frequency results. 

TFout_toavg = TFout_mixed;   %define here the time-frequency dataset to process.

twindow = [-0.2 1.0];
BaseLine = [twindow(1) -0.05];
TFout_avg = cell(size(TFout_toavg));

fdcfg = [];
fdcfg.keeptrials = 'no';    % find average across trials
fdcfg.channels = 'all';
fdcfg.frequency = 'all';
fdcfg.latency = twindow;
fdcfg.jackknife = 'yes';

selcfg = [];
selcfg.latency = BaseLine;
selcfg.avgovertime = 'yes';

selcfg2 = [];
selcfg2.latency = [0 twindow(2)];
selcfg2.avgovertime = 'no';

mcfg = [];
mcfg.parameter = 'powspctrm';
mcfg.operation = 'subtract';

for condcnt = 1:size(TFout_toavg,2)
    for sujcnt = 1:size(TFout_toavg,1)
        TFout_avg{sujcnt,condcnt} = ft_freqdescriptives(fdcfg, TFout_toavg{sujcnt,condcnt}); % Yields GA time-frequency data (Mixed) not baseline corrected
        blmean_curr = ft_selectdata(selcfg,TFout_avg{sujcnt,condcnt});
        poststim_curr = ft_selectdata(selcfg2,TFout_avg{sujcnt,condcnt});
        diffcurr = bsxfun(@minus,blmean_curr.powspctrm,TFout_avg{sujcnt,condcnt}.powspctrm);
        tvals_curr = diffcurr ./ TFout_avg{sujcnt,condcnt}.powspctrmsem;
        df_curr = size(tvals_curr,3)-1;             % the degrees of freedom (time points)
        pvals_curr = tcdf(tvals_curr,df_curr);
        zscore_curr = norminv(pvals_curr);
        TFout_avg{sujcnt,condcnt}.tvals = tvals_curr;
        TFout_avg{sujcnt,condcnt}.zscores = zscore_curr;
        TFout_avg{sujcnt,condcnt}.pvals = pvals_curr;
        TFout_avg{sujcnt,condcnt}.sigmap = [pvals_curr <= .01];
    end
end

clear condcnt sujcnt

%% Calculate the average z-score over all subjects.
% Based on calculation applied in http://www.jneurosci.org/content/28/34/8397#ref-30
elabels = {AllChans.chaninfo.labels};
selcfg = [];
selcfg.frequency = 'all';
selcfg.avgoverfreq = 'no';

GAzscore_foi = cell(size(TFout_avg));

for ccounter = 1:size(TFout_avg,2)
    for scounter = 1:size(TFout_avg,1)
        GAzscore_foi{scounter,ccounter} = ft_selectdata(selcfg,TFout_avg{scounter,ccounter});
    end
end

zgacfg = [];
zgacfg.keepindividual = 'no';
zgacfg.parameter = 'zscores';

GAzscore = cell(1,size(TFout_avg,2));

for condcntr = 1:size(TFout_avg,2)
    
    GAzscore{1,condcntr} = ft_freqgrandaverage(zgacfg,GAzscore_foi{:,condcntr});
    zscore_avg = sqrt(size(TFout_avg,1)).*GAzscore{1,condcntr}.zscores;
    GAzscore{1,condcntr}.zavg = zscore_avg;
    pvals_norm = normcdf(zscore_avg);
    pvals_norm2 = 2*normcdf(-abs(zscore_avg));   %two-tailed test
    pvals_lognorm = -log10(pvals_norm2);
    pvals_lognorm(isinf(pvals_lognorm)) = -log10(1);
    [GAzscore{1,condcntr}.pvals, GAzscore{1,condcntr}.sigmap] = fdr(pvals_norm2,.01);
    GAzscore{1,condcntr}.pvals = pvals_lognorm.*GAzscore{1,condcntr}.sigmap;
    
end

time = TFout_avg{1,1}.time;
freqs = GAzscore{1,1}.freq;
pres_type = 'map';      %or map

for condcntr2 = 1:size(TFout_avg,2)
    fhndl(condcntr2) = figure;
    set(gcf,'Name',TFcfg.Condnom1{1,condcntr2},'Color',[1 1 1]);
    for ecounter = 1:64
       subplot(8,8,ecounter)
       if strcmp(pres_type,'map')
           im(ecounter) = imagesc(time,freqs,squeeze(GAzscore{1,condcntr2}.pvals(ecounter,:,:)));
           set(gca,'YDir','normal') 
           title(elabels{1,ecounter})
       elseif strcmp(pres_type,'spectrum')
           plot(time,squeeze(GAzscore{1,condcntr2}.pvals(ecounter,:,:)));
            title(elabels{1,ecounter})
       end
    end
end

%% CARRY OUT BASELINE CORRECTION OF TIME-FREQUENCY DATA
% The following approaches are possible
% 1. Baseline correct at the level of single trials using subject-level
% mean baseline.
% 2. Baseline correct the mean of the subject-level GA time-frequency data
% using the subject-level mean baseline.
% 3. Baseline correct at the single trial level using the condition level
% mean baseline.
% 4. Baseline correct the mean of the subject-level GA time-frequency data
% using the condition-level mean baseline.
% 5. Baseline correct the single trial time-frequency data at the level of
% individual trials.

% EXTRACT THE MEAN SUBJECT-LEVEL REFERENCE INTERVAL (BASELINE)
% Replace the trial-level baseline by mean subject-level baseline.

blindx = [TFout_avg{1,1}.time >= BaseLine(1) & TFout_avg{1,1}.time <= BaseLine(2)];
blindx2 = [TFout_toavg{1,1}.time >= BaseLine(1) & TFout_toavg{1,1}.time <= BaseLine(2)];
TFout_toavg_mbl = TFout_toavg;

% Replace trial-level baseline by subject-level mean baseline.
for condcnt = 1:size(TFout_avg,2)
    fprintf('***condition number %d\n',condcnt);
    for sujcnt = 1:size(TFout_avg,1)
        fprintf('***subject %d\n',sujcnt)
        tnum = size(TFout_toavg{sujcnt,condcnt}.powspctrm,1);
        baseline_curr = repmat(TFout_avg{sujcnt,condcnt}.powspctrm(:,:,blindx),1,1,1,tnum);
        baseline_curr = permute(baseline_curr,[4 1 2 3]);
        TFout_toavg_mbl{sujcnt,condcnt}.powspctrm(:,:,:,blindx2) = baseline_curr;
    end
end

clear condcnt sujcnt

%% CORRECT THE TRIAL-LEVEL TIME-FREQUENCY DATA USING SUBJECT MEAN BASELINE
% AVERAGE ACROSS TRIALS 
TFout_trial_mbl = cell(size(TFout_toavg_mbl));   % Baseline corrected data at trial level
TFout_avg_mbl = cell(size(TFout_toavg_mbl));     % Subject-level mean tf data after baseline correction at single-trial level.
GAsub_avg_mbl = cell(1,size(TFout_toavg_mbl,2));
blc_type = 'zscore';

BLcfg = [];
BLcfg.baseline = BaseLine;
BLcfg.baselinetype = blc_type; 
BLcfg.parameter = 'powspctrm';

fdcfg = [];
fdcfg.keeptrials = 'no';    % find average across trials
fdcfg.channels = 'all';
fdcfg.frequency = 'all';
fdcfg.latency = twindow;

GAcfg = [];
GAcfg.keepindividual = 'yes';
GAcfg.foilim = 'all';
GAcfg.toilim = twindow;
GAcfg.channel = cellstr(TFout_avg{1,1}.label(1:72));
GAcfg.parameter = 'powspctrm';

TFGA_avg_mbl = cell(1,size(TFout_avg,2));
TFout_avg_nombl = cell(size(TFout_toavg_mbl)); 
TFGA_avg_nombl = cell(1,size(TFout_avg,2));

fprintf('***Carrying out single-trial baseline correction with subject-level mean baseline***');

for condcnt = 1:size(TFout_toavg,2)
    
    fprintf('***condition number %d\n',condcnt);
    for sujcnt = 1:size(TFout_toavg,1)
        fprintf('***subject %d\n',sujcnt) 
        TFout_trial_mbl{sujcnt,condcnt} = ft_freqbaseline(BLcfg,TFout_toavg_mbl{sujcnt,condcnt}); % Baseline corrected mixed time-frequeny data. 
        TFout_avg_mbl{sujcnt,condcnt} = ft_freqdescriptives(fdcfg, TFout_trial_mbl{sujcnt,condcnt}); 
        TFout_avg_nombl{sujcnt,condcnt} = ft_freqdescriptives(fdcfg, TFout_toavg_mbl{sujcnt,condcnt}); 
    end
    GAcfg.keepindividual = 'yes';
    TFGA_avg_mbl{1,condcnt} = ft_freqgrandaverage(GAcfg,TFout_avg_mbl{:,condcnt});
    TFGA_avg_nombl{1,condcnt} = ft_freqgrandaverage(GAcfg,TFout_avg_nombl{:,condcnt});  % Not baseline corrected
    GAcfg.keepindividual = 'no';
    GAsub_avg_mbl{1,condcnt} = ft_freqgrandaverage(GAcfg,TFout_avg_mbl{:,condcnt});
end

clear condcnt sujcnt

%% NORMALISE THE SUBJECT-LEVEL GA TIME-FREQUENCY DATA WITH SUBJECT-LEVEL MEAN BASELINE 

TFout_avg_bl = cell(size(TFout_avg));
TFGA_avg = cell(1,size(TFout_avg,2));
TFGA_avg_bl = cell(1,size(TFout_avg,2));

BLcfg = [];
BLcfg.baseline = BaseLine;
BLcfg.baselinetype = blc_type; 
BLcfg.parameter = 'powspctrm';

GAcfg = [];
GAcfg.keepindividual = 'yes';
GAcfg.foilim = 'all';
GAcfg.toilim = twindow;
GAcfg.channel = cellstr(TFout_avg{1,1}.label(1:72));
GAcfg.parameter = 'powspctrm';

fprintf('***Carrying out subject-level GA baseline correction with subject-level mean baseline***');

for condcnt = 1:size(TFout_avg,2)
    fprintf('***condition number %d\n',condcnt);
    for sujcnt = 1:size(TFout_avg,1)
        fprintf('***subject %d\n',sujcnt)
        TFout_avg_bl{sujcnt,condcnt} = ft_freqbaseline(BLcfg,TFout_avg{sujcnt,condcnt});
    end
    TFGA_avg{1,condcnt} = ft_freqgrandaverage(GAcfg,TFout_avg{:,condcnt});
    TFGA_avg_bl{1,condcnt} = ft_freqgrandaverage(GAcfg,TFout_avg_bl{:,condcnt});
end

%% PLOT THE RESULT OF THE JACKKNIFE METHOD 

cond_names = TFcfg.Condnom1;
eois = {'Fz' 'FCz' 'Cz' 'CPz'};
chindx = zeros(length(eois),1);
for cntr1 = 1:length(eois)
    i = cell2mat(cellfun(@strcmp,{elabels},{eois(1,cntr1)},'UniformOutput',false));
    chindx(cntr1) = find(i);
end

foi = [8:12];
selcfg1 = [];
selcfg1.frequency = [foi(1) foi(end)];

GAsub_avg_mbl_foi = cell(size(GAsub_avg_mbl));

for counter = 1:size(GAsub_avg_mbl,2)
 GAsub_avg_mbl_foi{1,counter} = ft_selectdata(selcfg1,GAsub_avg_mbl{1,counter});
end                                              

if ~strcmp(selcfg1.frequency,'all')
    findx = find([freqs>=foi(1) & freqs<=foi(end)]);
else
    findx = 1:length(freqs);
    foi = freqs;
end

for condcntr2 = 1:size(TFout_avg,2)
    
    fhndl(condcntr2) = figure;
    set(fhndl(condcntr2),'Name',TFcfg.Condnom1{1,condcntr2},'Color',[1 1 1]);
    for ecounter = 1:length(eois)
       subplot(2,2,ecounter)
       signind = squeeze(GAzscore{1,condcntr2}.sigmap(ecounter,findx,:));
       pow_curr = squeeze(GAsub_avg_mbl_foi{1,condcntr2}.powspctrm(chindx(ecounter),:,:));
       pow_curr(~istrue(signind)) = 0;
       im(ecounter) = imagesc(time,foi,pow_curr,[-4 4]);
       im(ecounter).UserData = {time,foi,pow_curr,elabels{1,chindx(ecounter)},cond_names{1,condcntr2}};
       colormap jet
       title(elabels{1,chindx(ecounter)})
       axhndl = gca;
       set(axhndl,'YDir','normal','NextPlot','add')
       xlabel('Time (secs)'); ylabel('Frequency (Hz)');
       set(fhndl(condcntr2),'CurrentAxes',axhndl);
       set(axhndl,'HitTest','on','SelectionHighlight','on','UserData',{time,freqs,pow_curr,elabels{1,chindx(ecounter)},cond_names{1,condcntr2}},...
           'ButtonDownFcn',{@plotsingleBoot_elec});
       
    end
end


%% CARRY OUT SINGLE TRIAL BASELINE CORRECTION OF TIME-FREQUENCY DATA

TFout_trialblc = cell(size(TFout_toavg));
TFout_trialblc_avg = cell(size(TFout_toavg));
TFGA_trialblc = cell(1,size(TFout_toavg,2));

BLcfg = [];
BLcfg.baseline = BaseLine;
BLcfg.baselinetype = blc_type; 
BLcfg.parameter = 'powspctrm';

fdcfg = [];
fdcfg.keeptrials = 'no';    % find average across trials
fdcfg.channels = 'all';
fdcfg.frequency = 'all';
fdcfg.latency = twindow;

GAcfg = [];
GAcfg.keepindividual = 'yes';
GAcfg.foilim = 'all';
GAcfg.toilim = twindow;
GAcfg.channel = cellstr(TFout_toavg{1,1}.label(1:64));
GAcfg.parameter = 'powspctrm';

fprintf('***Carrying out trial-level baseline correction ***');

for condcnt = 1:size(TFout_toavg,2)
    fprintf('***condition number %d\n',condcnt);
    for sujcnt = 1:size(TFout_toavg,1)
        fprintf('***subject %d\n',sujcnt)
        TFout_trialblc{sujcnt,condcnt} = ft_freqbaseline(BLcfg,TFout_toavg{sujcnt,condcnt});
        TFout_trialblc_avg{sujcnt,condcnt} = ft_freqdescriptives(fdcfg,TFout_trialblc{sujcnt,condcnt});
    end
    TFGA_trialblc{1,condcnt} = ft_freqgrandaverage(GAcfg,TFout_trialblc_avg{:,condcnt});
end

%% Baseline-correct the evoked data at the level of each subject. 
% Find the GA of the baseline corrected data but keep the individual
% subjects.

TFout_bl_evoked = cell(size(TFout_evoked));

blcfg = [];
blcfg.baseline = BaseLine;  %as TF analyNSIs do not chose 0 limit for baseline correction.
blcfg.baselinetype = blc_type;       % (data-mean_baseline)/mean_baseline
blcfg.parameter ='powspctrm';

for condcnt = 1:size(TFout_evoked,2)
    for sujcnt = 1:size(TFout_evoked,1)
        
        TFout_bl_evoked{sujcnt,condcnt} = ft_freqbaseline(blcfg, TFout_evoked{sujcnt,condcnt});
        
    end
end

TFGA_BL_ev = cell(1,size(TFout_evoked,2));
TFGA_ev = cell(1,size(TFout_evoked,2));

GAcfg=[];
GAcfg.keepindividual = 'yes';
GAcfg.foilim = 'all';
GAcfg.toilim = twindow;
GAcfg.channel = 'all';
GAcfg.parameter = 'powspctrm';

for condCnt=1:size(TFout_evoked,2)
    
    TFGA_BL_ev{1,condCnt}=ft_freqgrandaverage(GAcfg,TFout_bl_evoked{:,condCnt});  % Baseline corrected
    TFGA_ev{1,condCnt} = ft_freqgrandaverage(GAcfg,TFout_evoked{:,condCnt});      % Not baseline corrected
    
end

%% CALCULATE THE INDUCED ACTIVITY BY SUBTRACTING THE EVOKED DATA FROM THE MIXED DATA
% Here using the mixed time-frequency data in which the time-frequency
% decomposition was carried out at the single trial level and then averaged
% across trials before baseline correctio ==> subject-level mean baseline.

TFGA_BL_ind = cell(size(TFGA_BL_ev));
TFGA_noBL_ind = cell(size(TFGA_BL_ev));

for ccnt = 1:size(TFGA_BL_ev,2)
    
    TFGA_BL_ind{1,ccnt} = TFGA_BL_ev{1,ccnt};
    TFGA_noBL_ind{1,ccnt} = TFGA_ev{1,ccnt};
    powspct_diff = TFGA_avg_mbl{1,ccnt}.powspctrm - TFGA_BL_ev{1,ccnt}.powspctrm;
    powspct_diff_nobl = TFGA_avg_nombl{1,ccnt}.powspctrm - TFGA_ev{1,ccnt}.powspctrm;
    TFGA_BL_ind{1,ccnt}.powspctrm = [];
    TFGA_noBL_ind{1,ccnt}.powspctrm = [];
    TFGA_BL_ind{1,ccnt}.powspctrm = powspct_diff;
    TFGA_noBL_ind{1,ccnt}.powspctrm = powspct_diff_nobl;
end

%% VISUALISATION OF GRAND-AVERAGE TIME-FREQUENCY MAP

lcfg=[];
lcfg.layout= fullfile(filesep,'Users','bolger','Documents','MATLAB','fieldtrip-master','template','layout','biosemi64.lay');
lcfg.elec = TFout_toavg{1,1}.elec;
layout=ft_prepare_layout(lcfg,TFout_toavg{1,1});    %Prepare head layout
baseline = BaseLine;

figcfg=[];
figcfg.parameter='powspctrm';
figcfg.xlim = twindow;
figcfg.ylim='maxmin';
figcfg.zlim= [-8 8];
figcfg.masknans='yes';
figcfg.baseline= 'no';           %define the baseline for normalisation
figcfg.baselinetype ='';   %decibel converNSIon: specify power as change relative to baseline power. The correction is carried out on a NSIngle-trial baNSIs.
figcfg.hotkeys='yes';
figcfg.masktype='saturation';
figcfg.trials='all';
figcfg.box='yes';
figcfg.hotkeys='yes';
figcfg.colorbar='yes';
figcfg.showlabels='yes';
figcfg.interactive='yes';
figcfg.layout = layout;

for condCnt=1:length(condnames)
    figure; set(gcf,'Color',[1 1 1],'Name',condnames{condCnt,1});
    ft_multiplotTFR(figcfg,TFGA_avg_mbl{1,condCnt}) %plot all channels
end

%% Plot individual channel time-frequencies
chanois = {'Fz'; 'FCz'; 'Cz';'CPz'}; % 'Cz'; 'C4'; 'P3'; 'Pz'; 'P4'; 'PO3'; 'POz'; 'PO4'}; %;...
    %'CP1' 'CP3' 'CP5';'CP2' 'CP4' 'CP6';'P1' 'P3' 'P5';'P2' 'P4' 'P6'};

spcfg = [];
spcfg.parameter = 'powspctrm';
spcfg.baseline = 'no';
spcfg.baselinetype = 'zscore';
spcfg.xlim = [-0.2 1.0];
spcfg.ylim = 'maxmin';
spcfg.zlim = [-8 8];
spcfg.colormap = colormap('jet');
for condCnt=1:length(condnames)

    f1 = figure; 
    set(f1,'Color',[1 1 1],'Name',condnames{condCnt,1});
    for chancnt = 1:size(chanois,1)
        subplot(4,1,chancnt)
        spcfg.channel = {chanois{chancnt,:}};
        ft_singleplotTFR(spcfg, TFGA_BL_ind{1,condCnt}) %plot all channels
    end
end
%%
M = squeeze(statsout{1,1}.mask);
sigs = zeros(64,7);
M = M(:,1:size(time_curr,2));

for count = 1:64
    ind = find(M(count,:));
    if ~isempty(ind)
        sigs(count,ind) = time_curr(ind);
    end
    
end
%% 

topocfg = [];
topocfg.baseline = 'yes';
topocfg.parameter = 'powspctrm';
topocfg.baselinetype = 'zscore';
topocfg.xlim = [0.3 0.6];
topocfg.ylim = [8 12];
topocfg.zlim = [-10 10];
topocfg.layout = layout;
topocfg.colormap = 'jet';
topocfg.markers = 'off';
topocfg.highlight = 'off';

if datalen == 2
figure;
subplot(2,3,[1 2]); imagesc(time_curr,1:64,-log10(squeeze(statsout{1,1}.prob)))
subplot(2,3,4); ft_topoplotTFR(topocfg,TFGA_avg_nombl{1,2}); title(' Compatible (8-12Hz)');

subplot(2,3,5); ft_topoplotTFR(topocfg, TFGA_avg_nombl{1,1}); title('Incompatible (8-12Hz)');

topocfg.highlight = 'on';
cindx = find(sum(M,2));
topocfg.highlightchannel = cindx;

mcfg = [];
mcfg.parameter = 'powspctrm';
mcfg.operation = 'subtract';

Diff = ft_math(mcfg,TFGA_avg_nombl{1,2}, TFGA_avg_nombl{1,1});
subplot(2,3,6); ft_topoplotTFR(topocfg,Diff); title('Ironic - Literal (Non-sarcastic + Sarcastic): 4-7Hz');

elseif datalen ==4
    figure;
    subplot(2,3,[1 2]); imagesc(time_curr,1:64,-log10(squeeze(statsout{1,1}.prob)))
    
    m1cfg = [];
    m1cfg.parameter = 'powspctrm';
    m1cfg.operation = 'x1-x2';
    Diff1 = ft_math(m1cfg, DataIn{1,1}, DataIn{1,2});
    
    subplot(2,3,4); ft_topoplotTFR(topocfg,Diff1); title('Ironic (Non-sarcastic + Sarcastic): 11-13Hz');
    
    Diff2 = ft_math(m1cfg, DataIn{1,3}, DataIn{1,4});
    subplot(2,3,5); ft_topoplotTFR(topocfg, Diff2); title('Literal (Non-sarcastic + Sarcastic): 11-13Hz');
    
    topocfg.highlight = 'on';
    cindx = find(sum(M,2));
    topocfg.highlightchannel = cindx;
    
    m1cfg.operation = '(x1-x2)-(x3-x4)';
    Diff3 = ft_math(m1cfg,DataIn{1,1}, DataIn{1,2}, DataIn{1,3}, DataIn{1,4});
    
    subplot(2,3,6); ft_topoplotTFR(topocfg,Diff3); title('Ironic - Literal (Non-sarcastic): 11-13Hz');
    
    
end
%%
DataIn_sel = cell(1,datalen);
cfgsel = [];
cfgsel.channel = {'TP7' 'CP5' 'P9' 'P7'};  %}{AllChans.chaninfo(cindx).labels};
cfgsel.avgoverchan = 'yes';
cfgsel.latency = [0.3 0.6];
cfgsel.avgovertime = 'yes';
cfgsel.frequency = [11 13];
cfgsel.avgoverfreq = 'yes';

for dcnt = 1:datalen

    DataIn_sel{1,dcnt} = ft_selectdata(cfgsel, DataIn{1,dcnt});
    sem_curr(dcnt) = std(DataIn_sel{1,dcnt}.powspctrm)/sqrt(length(DataIn_sel{1,dcnt}.powspctrm));
    mean_curr(dcnt) = mean(DataIn_sel{1,dcnt}.powspctrm,1);
    sem_up(dcnt) = mean_curr(dcnt)+sem_curr(dcnt);
    sem_low(dcnt) = mean_curr(dcnt)-sem_curr(dcnt);
end


figure; 
b = bar(mean_curr,'FaceColor','flat');
b.CData(3:4) = 'g';
hold on
errorbar(1:datalen,mean_curr,sem_curr,'.k')

%% For individual electrode analysis

i= find(chan_alpha <= .05);
chanoi = {AllChans.chaninfo(i).labels};

figure; 
topoplot(-log10(chan_alpha),AllChans.chaninfo(1:64),'plotrad',0.5,'electrodes','on','plotchans',i,'maplimits',[0 3])

%% Carry out cluster-based permutation test to find differences between time-frequency results between two conditions

datalen = size(TFout_avg_mbl,2);
%datalen = size(DataIn,2);

if datalen ==2
    [data1, data2] = deal(TFGA_BL_ind{1,1}, TFGA_BL_ind{1,2});
    datanom1 = 'Incompatible';
    datanom2 = 'Compatible';
    dataIn = {data1, data2};
    datanomIn = {datanom1, datanom2};
elseif datalen ==4
    [data11, data12,data21,data22] = deal(DataIn{1,1}, DataIn{1,2},DataIn{1,3},DataIn{1,4});
    %[data11, data12,data21,data22] = deal(TFGA_avg_mbl{1,1}, TFGA_avg_mbl{1,2},TFGA_avg_mbl{1,3},TFGA_avg_mbl{1,4});
    datanom11 = 'NSI';
    datanom12 = 'NSL';
    datanom21 = 'SI';
    datanom22 = 'SL';
    dataIn = {data11,data12,data21,data22};
    datanomIn = {datanom11, datanom12, datanom21, datanom22};
end
          

% Need to include here option to investigate the mean over a band of
% frequencies.

cfg=[];
cfg.avgfreq = 'yes'; 
cfg.freqs = [11 13];    %[TFGA_avg_mbl{1,1}.freq(1:end)];
cfg.latency = [0.2 1];

if strcmp(cfg.avgfreq,'no')
    chan_oi = {AllChans.chaninfo(1:64).labels};
    cfg.DataIn = dataIn;
    cfg.mode = 'multielec'; %can also be indivelec multielec
    cfg.layout=layout;
    
    I = zeros(1,length(chan_oi));
    for cntr = 1:length(chan_oi)
        I(cntr) = find(strcmp({AllChans.chaninfo(1:64).labels},chan_oi{1,cntr}));
    end

    cfg.chanoi= {AllChans.chaninfo(I).labels}; 
    cfg.subnum = 25;
    cfg.CondNames = datanomIn;
    cfg.baseline = BaseLine;
    timeindx = diff(TFGA_avg_nombl{1,1}.time)*2;       %defined in seconds
    cfg.tsteps = [cfg.latency(1):timeindx:cfg.latency(2)-timeindx;cfg.latency(1)+timeindx:timeindx:cfg.latency(2)]';
   

elseif strcmp(cfg.avgfreq,'yes')

    selcfg.frequency   = cfg.freqs;
    selcfg.avgoverfreq = 'yes';
    selcfg.latency = cfg.latency;
    
    Data_in = cell(1,size(dataIn,2));
    
    for x = 1:size(dataIn,2)
        Data_in{1,x} = ft_selectdata(selcfg, dataIn{1,x});    
    end
   
    cfg.freqs = Data_in{1,1}.freq;
    cfg.DataIn = dataIn;              % overwrite the input data by providing just the mean over the frequencies specified.
    cfg.mode='multielec';              % can also be indivelec multielec
    cfg.layout=layout;
    cfg.chanoi = {AllChans.chaninfo(1:64).labels};
    cfg.subnum = 25;
    cfg.CondNames = datanomIn;
    cfg.baseline = BaseLine;
    timeindx = diff(TFGA_avg_nombl{1,1}.time)*6;       %defined in seconds
    cfg.tsteps = [cfg.latency(1):timeindx:cfg.latency(2)-timeindx;cfg.latency(1)+timeindx:timeindx:cfg.latency(2)]';
    

end

Out_stat = CREx_TFstats_cluster(cfg);

%% 

lcfg=[];
lcfg.layout='/Users/bolger/Documents/MATLAB/fieldtrip-master/template/layout/biosemi64.lay';
layout=ft_prepare_layout(lcfg,TFout_mixed{1,1});


CREx_TFactive_bl_stat(TFGA_noBL_ind{1,2},TFGA_BL_ind{1,2},layout,AllChans);

%% NEED TO EXTRACT FREQUENCY BAND DATA AND WRITE TO AN EXCEL FILE FOR FURTHER ANALYSIS - Projet Ironie
% ROI1=left-anterior-ventral (AF7, F7, FT7, F5, FC5, C5, T7), ROI2=left-anterior-dorsal (AF3, F3, FC3, F1, FC1, C1, C3), 
% ROI3=left-posterior-ventral (TP7, CP5, P5, P7, P9, PO7), and ROI4=left-posterior-dorsal (CP3, CP1, P3, PO3, P1, O1), 
% ROI5=midline-anterior (AFz, Fz, FCz, Cz), ROI6=midline-posterior (CPz, Pz, POz,Oz), 
% ROI7=right-anterior-ventral (AF8, F8, FT8, F6, FC6, C6, T8), ROI8=right-anterior-dorsal (AF4, F4, FC4, F2, FC2, C2, C4), 
% ROI9=right-posterior-ventral (TP8, CP6, P6, P8, P10, PO8), ROI10=right-posterior-dorsal (CP4, CP2, P4, PO4, P2, O2). 
pname = fullfile(filesep,'Users','bolger','Documents','work','Projects','projet-Ironie','Ironie-data',filesep);
freqband_oi = [55 80];      %Specify the frequency band of interest.
timewind_oi = [0.5 0.8];      %Specify the time-window of interest.
elec_groupsnoms = {'left-anterior-ventral','left-anterior-dorsal','left-posterior-ventral','left-posterior-dorsal',...
    'midline-anterior','mindline-posterior','right-anterior-ventral','righte-anterior-dorsal','right-posterior-ventral',...
    'right-posterior-dorsal'};
sujnum = 22;

elec_grouplabs = {{'AF7', 'F7', 'FT7', 'F5', 'FC5', 'C5', 'T7'},{'AF3', 'F3', 'FC3', 'F1', 'FC1', 'C1','C3'},{'TP7', 'CP5', 'P5', 'P7', 'P9', 'PO7'},...
    {'CP3', 'CP1', 'P3', 'PO3', 'P1', 'O1'},{'AFz', 'Fz', 'FCz', 'Cz'},{'CPz', 'Pz', 'POz','Oz'},{'AF8', 'F8', 'FT8', 'F6', 'FC6', 'C6','T8'},...
    {'AF4', 'F4', 'FC4', 'F2', 'FC2', 'C2', 'C4'},{'TP8', 'CP6', 'P6', 'P8', 'P10', 'PO8'},{'CP4', 'CP2', 'P4', 'PO4', 'P2', 'O2'}};

conds = {'NSI','NSL','SI','SL'};

seloutcfg = [];
seloutcfg.frequency = freqband_oi;
seloutcfg.avgoverfreq = 'yes';


for ecnt = 1:size(elec_grouplabs,2)
    
    powspect_curr = cell(1,size(TFGA_avg_mbl,2));
    roinames_curr = elec_groupsnoms{1,ecnt};
    elecnames_curr = elec_grouplabs{1,ecnt};
    
    for condcnt = 1:size(TFGA_avg_mbl,2)
        
        seloutcfg.channel = elec_grouplabs{1,ecnt};
        seloutcfg.avgoverchan = 'no';
        seloutcfg.latency = timewind_oi;
        seloutcfg.avgovertime = 'yes';
        
        DataCurr = ft_selectdata(seloutcfg,TFGA_avg_mbl{1,condcnt});
        powspect_curr{1,condcnt} = DataCurr.powspctrm;

    end
    % Write to excel file here...
    fnom = [num2str(freqband_oi(1)),'-',num2str(freqband_oi(2)),'Hz-',elec_groupsnoms{1,ecnt},'.txt'];
    fnomfull = [pname,fnom];
    sheet_curr = elec_groupsnoms{1,ecnt};
    powspect_mat = cell2mat(powspect_curr); 
    csvwrite(fnomfull,powspect_mat,3,1);
    
end

%In every Excel page there will be the mean over the frequency band of
%interest and over the time window of interest for each electrode
%comprising each of the 10 ROI for each of the 4 conditions. 











