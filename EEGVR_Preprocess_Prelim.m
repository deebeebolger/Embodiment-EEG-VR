%% Function to carry out the preliminary pre-processing steps of EEG-VR data.
% Date: May 2018                                Programmed by: D. Bolger
%****************************************************************************

%% DEFINE PATHS FOR SUBJECT LEVEL FILES, CHANNEL LABEL FILE AND FILE CONTAINING ELECTRODE COORDINATES.
% These paths need to be changed by a new user. 

% Path to base directory containing the subject-level data.
DIR_main = fullfile('Users','bolger','Documents','work','Projects','Project-EEG-VR','ExpData',filesep);

% Path to file with channel labels.
chandir = fullfile('Users','bolger','Documents','work','Projects','Project-EEG-VR','ExpData','Chaninfo.mat');

% Path to file containg electrode coordinates for 64-electrode 10-20 system.
chlocpath = fullfile('Users','bolger','Documents','MATLAB','eeglab13_6_5b','plugins','dipfit2.3',...
                'standard_BESA','standard-10-5-cap385.elp');  

% Path to xls file containing word-frequency and uniqueness point data. 
xls_wordfreq = fullfile('Users','bolger','Documents','work','Projects','Project-EEG-VR','EEG-VR_verb_frequencies.xlsx');

% Path to the wav files (auditory stimuli). 
wav_path = fullfile('Users','bolger','Documents','work','Projects','Project-EEG-VR','VR-Embodiment wav files','soundfiles-fr-final',filesep);

% Path to the parameters *.mat file. 
params_path = fullfile('Users','bolger','Documents','work','Projects','Project-EEG-VR','ExpData','EEG-VR_parameters.txt');
%% OPEN DIALOGUE BOX TO SPECIFY SUBJECT NUMBER AND LIST TO PROCESS.
prompt1={'Subject Number/s:' 'List/s:'};
dlg_title = 'Specify subject and list:';
num_lignes = [10;10];
prompt1_ans = inputdlg(prompt1,dlg_title,num_lignes);

subs=cellstr(prompt1_ans{1,1});
sujindx=cellfun(@str2double,subs);
listcurr = cellstr(prompt1_ans{2,1});

%% OPEN EEGLAB SESSION
% Opens an EEGLAB session and presents main EEGLAB GUI. 

addpath(genpath(''));
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


for counter = 1:length(sujindx)  
    %% GENERATE THE CURRENT SUBJECT-LEVEL FILE-PATH
    % Tries to generate a file path string to matlab the real file-path for current subject.

    if sujindx < 10
        subjnom = strcat('s0',num2str(sujindx(counter)),listcurr{counter,1});
    else
        subjnom = strcat('s',num2str(sujindx(counter)),listcurr{counter,1});
    end

    Dirbase = fullfile(DIR_main,subjnom,filesep);     % Define the current subject-level filepath.

%% LOAD IN THE *.bdf FILE OF THE CURRENT SUBJECT
% Loads in the *.bdf file, ensuring that 74 channels are included so that the 
% ERGO1 and ERGO2 data are loaded.

    allfiles= dir(Dirbase);
    fileIndex = find(~[allfiles.isdir]);
    filenum=dir(strcat(Dirbase,'*.bdf'));                      %find all the *.bdf files in the current folder
    filenom={filenum.name};


    fullDir = strcat(Dirbase,filenom{1,1});
    fnom = subjnom;

    % The following three lines is added to resolve a bug occurring when
    % opening the *.bdf file.
    x = fileparts( which('sopen') );
    rmpath(x);
    addpath(x,'-begin');

    % Opening up *.bdf file and saving as a *.set file.
    EEG = pop_biosig(fullDir, 'channels',[1:74], 'ref', [] ,'refoptions',{'keepref' 'off'} );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom),'gui','off'); % Create a new dataset for the current raw datafile
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom),'filepath',Dirbase);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw


%% CALL OF FUNCTION TO EXTRACT TRIAL DATA OUTPUT BY THE STIMULATION SYSTEM.
% The trial level files output by unity should be in a "StimData/Trials/" folder saved in the current subject's folder.
% A stimdata.txt file will be saved in the subject-level directory.
% However, if this stimdata.txt file already exists for the subject,
% this step is skipped. 

    stimdir = fullfile(Dirbase,'StimData',filesep);
    xtest = dir(stimdir);
    filesnom = {xtest.name};

    stim_test = strcat(subjnom,'_stimdata.txt');
    if sum(ismember(filesnom,stim_test))==0
        disp('-------Create stimdata file--------');
        stimfile = EEGVR_extract_trial_data(subjnom,Dirbase); % Call of function to assemble trial-level Unity output files. 
    else
        disp('-------Stimdata file already present--------');
    end

%% EXTRACT THE ERGO DATA (AUDITORY STIMULI AND PHOTODIODE SIGNALS) 
% the photodiode and auditory stimuli channels are extracted from the EEG and saved to a separate
% field of the EEG structure, "EEG.dataERGO".
% It, therefore, erases the channels 73 and 74 and saves the changes.
% A figure of the auditory stimuli and the photodiode signals over is
% presented. 

    EEG.dataERGO = EEG.data(73:74,:);

    EEG = pop_select( EEG,'nochannel',{'Erg1' 'Erg2'});
    EEG = eeg_checkset( EEG );

    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom),'filepath',Dirbase);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw

    f1 = figure;
    subplot(2,1,1); plot(EEG.times./1000,EEG.dataERGO(1,:))
    title(strcat('Photodiode Signals:',subjnom));
    xlabel('Time (seconds)');
    subplot(2,1,2); plot(EEG.times./1000,EEG.dataERGO(2,:));
    title(strcat('Auditory Stimuli:',subjnom));
    xlabel('Time (seconds)');


%% ADD CHANNEL INFORMATION TO THE CURRENT DATASET.
% Channel coordinates and labels are added to the current dataset and
% the dataset is saved to the current subject-level directory. 
% The Chaninfo.mat file is loaded as it contains the electrode labels. 
% From EEGLAB plugins, the file, "standard-10-5-cap385.elp" is loaded
% as this contains the correct coordinates for the 10-20 system used here. 

    chaninfo = load(chandir,'Chaninfo');
    chans = chaninfo.Chaninfo;

    for cnt = 1:length(chans)
        EEG.chanlocs(cnt).labels = chans{1,cnt};
    end

    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );

    % !!!! For a new user, need to change this file path.
    
    EEG=pop_chanedit(EEG, 'lookup',chlocpath);                % Load channel path information
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom),'filepath',Dirbase);
    eeglab redraw

%% NEED TO ADD EVENT INFORMATION BASED ON DATA RECORDED ONTO ERGO CHANNELS (AUDITORY STIMUI AND PHOTODIODE SIGNAL).
% Trig_onsets_detect() function finds photodiode onset times.
% EEGVR_detect_audonsets() function finds the auditory stimuli onset times.

    trialtrigs = Trig_onsets_detect(EEG.dataERGO,EEG.times,EEG.srate,Dirbase);  % Call function to identify photodiode trigger onsets.
    audonsets = EEGVR_detect_audonsets(EEG.dataERGO(2,:),EEG.times,EEG.srate,Dirbase,subjnom,trialtrigs,wav_path);  % Call function to identify auditory stimuli onsets

    % Merging the trialtrigs and audonsets variables.
    T = trialtrigs(~cellfun('isempty',trialtrigs));
    Trialtrigs = reshape(T,[length(T)/size(trialtrigs,2),size(trialtrigs,2)]);
    ttrigs = cell(size(Trialtrigs,1),size(Trialtrigs,2));

    for cnt = 1:length(Trialtrigs)

        ttrigs{cnt,1} = Trialtrigs{cnt,1};
        ttrigs{cnt,2} = Trialtrigs{cnt,2};
        ttrigs{cnt,3} = Trialtrigs{cnt,3};
        ttrigs{cnt,4} = Trialtrigs{cnt,4};
    end

    % Sort latency data for auditory stimulus and photodiode signal onsets.
    lat_all = cat(1,audonsets{:,4},ttrigs{:,4});
    time_all = cat(1,audonsets{:,1},ttrigs{:,2});
    type_all = cat(1,{audonsets{:,2}}',{ttrigs{:,1}}');
    [latencies_all,indx] = sort(lat_all,'ascend');
    lat_all = lat_all(indx);
    time_all = time_all(indx);
    type_all = {type_all{indx}}';

% Integrate the onset times into the event field of the EEG structure.
%Create the EEG.event field from the "audonsets" and "Trigtrials"
% variables.
    EEG.event = [];
    EEG.urevent = [];

    for cnt1 = 1:length(type_all)
        EEG.event(cnt1).type = type_all{cnt1,1};
        EEG.event(cnt1).latency = lat_all(cnt1,1);
        EEG.event(cnt1).urevent = cnt1;
        EEG.urevent(cnt1).type = type_all{cnt1,1};
        EEG.urevent(cnt1).latency = lat_all(cnt1,1);
    end

    % Save the dataset with onset information added to the current
    % subject-level directory. 
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom),'filepath',Dirbase);
    eeglab redraw

%% ASSIGN TRIAL-TYPES (GO OR NOGO) AND INFORMATION REGARDING "GOOD" AND "BAD" TRIALS.
% This facilitates epoching only go or no-go trials separately and the
% automatic rejection of bad trials after epoching. 
% Note that the word "bad" is added to auditory stimuli of incorrect
% trials. 

    [EEG,currconds] = EEGVR_trialtype_assign(EEG,subjnom,Dirbase);

    if ~exist('fnom')      % Sometimes the dataset does not have a setname field defined. 
        fnom = EEG.setname;
    end

    % Save the dataset with this trial information added to the current
    % subject-level directory. 
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom),'filepath',Dirbase);
    eeglab redraw

%% PREPARE INFORMATION TEXT-FILE FOR THE CURRENT SUBJECT. 

    fname=strcat(fnom,'-info.txt');
    fdir=strcat(Dirbase,fname);
    fid=fopen(fdir,'w');
    fprintf(fid,['---------',fnom,'----------\n\n']);

    % Adds information regarding the total number of trials and the total number of "go","nogo" and
    % "bad" trials to this text file. 
    nogo_num = length(find(strcmp(currconds{1,5},'NOGO')));     % The number of "nogo" trials. 
    go_num = length(currconds{1,5}) - nogo_num;                 % The number of "go" trials.
    bad_num = length(find(strcmp(currconds{1,7},'badtrial')));  % The number of bad trials. 
    verbs = currconds{1,2};
    badwords = verbs(strcmp(currconds{1,7},'badtrial'));        % The verbs corresponding to the bad trials. 
    badtrial_indx = find(strcmp(currconds{1,7},'badtrial'));
    
    fprintf(fid,'The total number of trials is: %d\n',length(currconds{1,5}));
    fprintf(fid,'The total number of GO trials is: %d\n',go_num);
    fprintf(fid,'The total number of NOGO trials is: %d\n',nogo_num);
    fprintf(fid,'The number of incorrect trials: %d\n',bad_num);
    fprintf(fid,'The verbs corresponding to incorrect trials:');
    for i = 1:bad_num
        fprintf(fid,'%s\t',badwords{i,1});
    end
    fprintf(fid,'\nTrial indices corresponding to incorrect trials:');
    for i1 = 1:bad_num
        fprintf(fid,'%d\t',badtrial_indx(i1));
    end
    fprintf(fid,'\n---------------------------------------\n'); 

%% OPEN THE PARAMETERS FILE WITH PRE-PROCESSING PARAMETERS.
% Load the pre-processing parameters defined in the parameters *.txt
% file into the current workspace. 

    % !!! Need to change this file path for a new user.
    fid2=fopen(params_path);      % il faut changer le chemin
    mydata=textscan(fid2,'%s %s');

    for i = 1:length(mydata{1,1})                     % generate a parameters structure from the parameters text file
        Params.(genvarname(mydata{1,1}{i}))=mydata{1,2}(i);
    end

    f_low = str2double(Params.fc_low);
    f_hi= str2double(Params.fc_hi);
    SR = str2double(Params.srate); 
    SR_orig = EEG.srate; 

%% RESAMPLE THE CONTINUOUS DATA USING THE PARAMETER DEFINED IN THE PARAMETERS TEXT FILE. 
% Resamples using the EEG resampling function. 
% If the user has the Matlab signal processing toolbox, it uses the
% Matlab resample() function. 
% Write information regarding the resampling to the subject-level text
% file. 

    fprintf(fid,'\nDownsampled from %fHz to %fHz\n',SR_orig,SR);
    display('***********************************Resampling to 512Hz*******************************************')
    fnom_rs = strcat(fnom,'-rs');

    EEG = pop_resample(EEG, SR);   %resample the data at sampling rate defined, sr.
    EEG =eeg_checkset(EEG);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_rs),'gui','off'); % current set = xx;
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom_rs),'filepath',Dirbase);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw


%% APPLY BAND-PASS FILTER BETWEEN THE LOWER AND UPPER LIMITS SPECIFIED IN PARAMETERS FILE.
% It applies a FIR windowed sinc filter using a blackman filter.
% The filter frequency response is plotted. 
% The details of the filtering are written to subject information txt file.

    display('*********************************Bandpass filtering using a FIR windowed sinc filter***********************************')
    [M, dev]=pop_firwsord('blackman',SR, 2);
    [EEG,com,b]=pop_firws(EEG,'fcutoff',[f_low f_hi],'forder',M,'ftype','bandpass','wtype','blackman');
    fvtool(b);                                      % Visualise the filter characteristics
    fnom_filt=strcat(fnom_rs,'-filt');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_filt),'gui','off');   %save the resampled data as a newdata set.
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom_filt),'filepath',Dirbase);
    eeglab redraw

    fprintf(fid,'Band-pass filtered %f2 - %dHz with %d-order fir windowed sinc filter (blackman).\n',f_low,f_hi,M);

%% RE-REFERENCE THE DATA TO THE ELECTRODES SPECIFIED IN THE PARAMETERS FILE.
% The channels used for referencing are generally EXG1 and EXG2,
% channels 65 and 66, respectively. 
% The details of the re-referencing are written to the information text file. 

    R=num2str(Params.references{1,1});
    if length(R)==2
        refs=str2double(R);
    elseif length(R)==4
        refs =[str2double(R(1:2)) str2double(R(3:4))];
    end

    display('***********************Rereference to Defined Channel:  does zero potential exist?*****************************')
    EEG =pop_reref(EEG, refs, 'method','standard','keepref','on');
    fnom_ref = strcat(fnom_filt,'-rref');
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_ref),'gui','off');   %save the resampled data as a newdata set.
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom_ref),'filepath',Dirbase);
    EEG = eeg_checkset( EEG );
    eeglab redraw

    here = CURRENTSET;   % Mark the current set. 

    fprintf(fid,'Rereferenced using channels %s and %s.\n',EEG.chanlocs(refs(1)).labels,EEG.chanlocs(refs(2)).labels);
    fprintf(fid,'-------------------------------------------------------------------------------\n');

%% EXTRACT THE FIRST 30 SECONDS OF REST-STATE EEG FROM THE RESAMPLED, FILTERED AND RE-REFERENCED CONTINUOUS DATA.
% The user needs to specify the upper and lower time limits for the
% resting state data to avoid include data with movement artifacts etc.
% The resting-state data is saved to a separate *.set file in the
% current subject-level folder. 

    disp('----------------Extracting Resting State Interval-----------------------------');
    timex_lims = [10 40];  % Time limits in seconds. 
    EEG = pop_select( EEG,'time',timex_lims );
    fnom_rest = strcat(fnom_filt,'-resting');
    EEG.setname = fnom_rest;
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_rest),'gui','off');   %save the resampled data as a newdata set.
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom_rest),'filepath',Dirbase);
    EEG = eeg_checkset( EEG );
    eeglab redraw; 

%% DETECT POSSIBLE BAD ELECTRODES AUTOMATICALLY VIA EEGLAB KURTOSIS MEASURE.
% Those electrodes with a kurtosis value >5 (z-score) are marked as
% bad.
% Bad electrodes electrodes detected with the measure are written to
% the subject information text file.

    % Retrieve the dataset before the resting state. 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',here,'study',0);
    EEG = eeg_checkset( EEG );
    eeglab redraw

    chans = {EEG.chanlocs.labels};
    [EEG, eindx, measure,~] = pop_rejchan(EEG, 'elec',[1:64] ,'threshold',5,'norm','on','measure','kurt');
    EEG.reject.rejkurtE = eindx;                          %indices of suggested electrodes to reject according to kurtosis.

    if ~isempty(eindx)
        for cntr=1:length(eindx)
            if cntr==1
                fprintf(fid,'Bad electrodes according to kurtosis (threshold z-score 5):  %s  ',chans{eindx(cntr)});
            elseif cntr>1 && cntr<length(eindx)
                fprintf(fid,' %s  ',chans{eindx(cntr)});
            elseif cntr==length(eindx)
                fprintf(fid,' %s \n ',chans{eindx(cntr)});
            end

        end
    else
        fprintf(fid,'No bad electrodes marked according to kurtosis (threshold z-score 5)\n');
    end
    %% RUN THE PREP PIPELINE FUNCTION TO AUTOMATICALLY DETECT NOISY CHANNELS, findNoisyChannels() 
    % Before running this function will need to take out the EXG channels that do not have X Y Z coordinates. 
    % This dataset is only used for the noisy channel detection script. 
    % This PREP function applies 4 different measures:
    % 1. Robust standard deviation (unusually high or low amplitude)
    % 2. Signal-to-Noise Ratio (SNR) (Christian Kothes, clean_channels()
    % function).
    % 3. Global correlation criteria (Nima Bigdely-Shamlo).
    % 4. RANSAC correlation (but may not always be performed).

    EEG = pop_select( EEG,'nochannel',{'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8'});
    EEG.setname = strcat(EEG.setname,'-remove-EXG');
    EEG = eeg_checkset( EEG );

    noisyOut = findNoisyChannels(EEG);   % Call of PREP pipeline function. 
    
    badchan_indx = noisyOut.noisyChannels.all;
    badchans = {chans{[badchan_indx]}};

    fprintf(fid,'Bad channels according to PREP pipeline noisy channel detector:\n');
    for i2 = 1:length(badchans)
        fprintf(fid,'%s\t',badchans{1,i2}); 
    end

    % Retrieve the original correct dataset before suppression of EXG channels. 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',here,'study',0);
    EEG = eeg_checkset( EEG );
    eeglab redraw


%% EPOCH THE CONTINUOUS DATA BUT INCLUDING ALL CONDITIONS.
% Include all conditions (go and nogo) and all words.
% Epoching is carried out with the verb onset as the T0. 
    
    dirsave = EEG.filepath; 
    Enom = strcat(EEG.setname,'-allconds');
    Ecurr = {EEG.event.type};
    toexcl = {'right' 'left' 'go' 'nogo'};
    condindx = ~ismember(Ecurr,toexcl);
    
    EEG = pop_epoch( EEG, Ecurr(condindx), [str2double(Params.wind_low{1,1}) str2double(Params.wind_hi{1,1})], 'newname', char(Enom), 'epochinfo', 'yes');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Enom),'gui','off');
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(Enom),'filepath',dirsave);
    EEG = eeg_checkset( EEG );
    eeglab redraw; 
    
    % CARRY OUT BASELINE CORRECTION 
    disp('--------------------Baseline correction-----------------------------------');
    Enom_bl = strcat(Enom,'-bl');
    EEG = pop_rmbase( EEG, [Params.wind_low{1,1} 0]);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Enom_bl),'gui','off');
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(Enom_bl),'filepath',dirsave);
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
%% CALL OF FUNCTION TO ADD COVARIATE DATA (WORD FREQUENCY AND UP-PHON) TO THE EVENTS. 
% It is best if the covariate data is added to the epoched data, before
% cleaning. 

    [EEG] = EEGVR_addcovariates(EEG,Ecurr(condindx),xls_wordfreq);

    EEG = pop_saveset( EEG, 'filename',char(Enom_bl),'filepath',dirsave);
    EEG = eeg_checkset( EEG );
    eeglab redraw;

%% REMOVE INCORRECT TRIALS AUTOMATICALLY.
% The variable : batrial_indx(). These trials are then removed
% automatically.
    
    disp(horzcat('-----------Remove the ',num2str(length(badtrial_indx)),' incorrect trials----------------------'));
    EEG.badtrialindx = badtrial_indx;
    EEG = pop_select(EEG, 'notrial', badtrial_indx);
    EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',dirsave);
    EEG = eeg_checkset( EEG );
    eeglab redraw;

    
%% AUTOMATIC BAD TRIAL DETECTION.
% Applies a threshold of 75mV and a maximum of 10% rejection.
        
    fprintf(fid,'\n*******************Condition %s:\n %d trials before cleaning\n','Go-NoGo',size(EEG.data,3));
    [EEG,remepochs] = pop_autorej(EEG, 'nogui','on','threshold',75,'eegplot','on','maxrej',10);
    numrej = find(EEG.reject.rejauto);
    EEG.trialrejauto = remepochs;
    
    EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',dirsave);
    EEG = eeg_checkset( EEG );
    eeglab redraw;
    
    fprintf(fid,'Number of trials marked for rejection (>75mV): %d\n', length(numrej));
    fprintf(fid,'Indx of trials marked for rejection (>75mV):\t');
    for cnt = 1:length(remepochs)
        fprintf(fid,'%d\t',remepochs(cnt))
    end
    fprintf(fid,'\n');
    fclose(fid);

%% CALCULATE THE SPECTRUM OF EACH ELECTRODE USING MULTI-TAPER (dpss). 
% Plot the mean spectrum of each electrode over all trials.
% The figure is interactive: if you select an electrode, its label will be
% prented in the command window. 
% Looking at the spectra of the electrodes may detect noisy electrodes. 

   CREx_SpectCalc_multitap(EEG,chans);

%% RUN FUNCTION TO LOCATE BAD EPOCHS AND CHANNELS VISUALLY.
    
   EpochChan_dlg(EEG); 

end % end of sujindx counter loop






