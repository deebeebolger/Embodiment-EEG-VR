function FTIn_rs=CREx_EEGFT_configure(EEGIn,fpath,fnom, rsopt,targevnt)
% Date: 07-09-2017   Programmed by: D. Bolger
% Function to configure EEGLAB data for use in FieldTrip toolbox
% Input data:
% EEGIn : Input EEGLAB structure
% DirIn : Directory of EEGIn data
% padlims: the time at the start and end of each trial to which to extend
% the data.
% targevnt - specify 'code' of a target event that you wish to isolate
%******************************************************
fulldir = fullfile(fpath,fnom);

assignin('base','EEGIn',EEGIn)
assignin('base','fulldir',fulldir);

if isempty(EEGIn.filepath) || isempty(EEGIn.filename)
    infile = fullfile(fulldir);
    hdr=ft_read_header(infile);
    events = ft_read_event(infile);
else
    infile = fullfile(EEGIn.filepath,EEGIn.filename);   %path of the *.set file
    hdr=ft_read_header(infile);
end 

dataIn = ft_read_data(infile,'header',hdr);
ftype = 'eeglab_set';

% Define trial information
if ischar([EEGIn.event.type])==1
    disp('Triggers are in string form')
    eventValues=unique({EEGIn.event.type});
else
    eventValues=unique([EEGIn.event.type]);
end

if ~isempty(targevnt)
    Bid = contains(eventValues,targevnt);
    TargetE = eventValues(1,Bid);
else
    TargetE = eventValues;
end

assignin('base','TargetE',TargetE);

prepost = [EEGIn.times(1)/1000 EEGIn.times(end)/1000];
cfg = [];
cfg.data=infile;
cfg.dataset = infile;
cfg.headerfile = infile;
cfg.dataformat = ftype;
cfg.headerformat = ftype;
cfg.continuous = 'no';                     %segmented data
cfg.trialdef.prestim = abs(prepost(1));    %note that it is the abs value...
cfg.trialdef.poststim = prepost(2);
cfg.trialdef.eventtype = 'trigger';
cfg.trialdef.eventvalue = TargetE ;
cfg.trialdef.ntrials = size(EEGIn.data,3);
cfg.channel='all';

tdef = ft_definetrial(cfg);
tdef.channel='all';
FTIn = ft_preprocessing(tdef);
FTIn_rs=[];

if strcmp(rsopt,'yes')
    
    cfg=[];
    cfg.resamplefs=256;
    cfg.detrend='no';
    cfg.demean='no';
    cfg.feedback='text';
    cfg.trials='all';
    
    FTIn_rs=ft_resampledata(cfg,FTIn);
    
else
    FTIn_rs=FTIn;
end


end % end of function