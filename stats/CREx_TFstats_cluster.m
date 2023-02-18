function statsout=CREx_TFstats_cluster(cfg)
% Date: 03-06-2016                       Programmed: D. Bolger
% The DataIn variable should be a 1x2 cell arry containing the data of the
% two experimental conditions.
% Mode : 'indivelec' for cluster-based permutation test on individual
%             specified electrodes
%            'multielec' for cluster-based permutation test on
%            channel-time-frequency data.
% layout : layout data to be used to prepare the neighbour data.
% (structure)
% chans:  labels of the channels over which to carry out analysis (cell
% array of strings)
% subnum: the number of subjects
% foi : frequency/frequencies of interest
% bline: baseline limits in seconds: [min max]
%***************************************************************************************

%PREPARE DATA
if size(cfg.DataIn,2) == 2
    [CondA, CondB] = deal(cfg.DataIn{1,1},cfg.DataIn{1,2});
    tindx=find([CondA.time>=cfg.latency(1) & CondA.time<=cfg.latency(2)]);
    time_curr = CondA.time(tindx);
elseif size(cfg.DataIn,2) ==4
    [CondA1, CondA2, CondB1, CondB2] = deal(cfg.DataIn{1,1},cfg.DataIn{1,2},cfg.DataIn{1,3},cfg.DataIn{1,4});
    tindx=find([CondA1.time>=cfg.latency(1) & CondA1.time<=cfg.latency(2)]);
    time_curr = CondA1.time(tindx);
    assignin('base','CondA1',CondA1)
    assignin('base','CondA2',CondA2);
    assignin('base','CondB1',CondB1);
    assignin('base','CondB2',CondB2);
end

Mode = cfg.mode;

layout = cfg.layout;
chans = cfg.chanoi;
foi = cfg.freqs ;
subnum=cfg.subnum;
% bline = cfg.baseline;
Condnom = cfg.CondNames;
freqavg = cfg.avgfreq;
latency = cfg.latency;
time_steps = cfg.tsteps;

% Find the index of the total trial time and the time index for current
% latency.


assignin('base','time_curr',time_curr);
assignin('base','tindx',tindx);


%% CALCULATE THE RAW EFFET

% Find the indices of the frequencies of interest
% Compare CondA.freq and foi.
%length of channels here is generally all or 64 channels.

ind=zeros(length(foi),1);
raweffect=cell(length(foi),1);
statsout=[];

if size(cfg.DataIn,2) ==2
    meanCondA=squeeze(nanmean(CondA.powspctrm(:,1:length(chans),:,tindx),1));  %dimensions 64 x numfreq x timepoints
    meanCondB=squeeze(nanmean(CondB.powspctrm(:,1:length(chans),:,tindx),1));
    assignin('base','meanCondA',meanCondA);
    assignin('base','meanCondB',meanCondB);
elseif size(cfg.DataIn,2) == 4
    meanCondA1 = squeeze(nanmean(CondA1.powspctrm(:,1:length(chans),:,tindx),1));  %dimensions 64 x numfreq x timepoints
    meanCondA2 = squeeze(nanmean(CondA2.powspctrm(:,1:length(chans),:,tindx),1));
    meanCondB1 = squeeze(nanmean(CondB1.powspctrm(:,1:length(chans),:,tindx),1));  %dimensions 64 x numfreq x timepoints
    meanCondB2 = squeeze(nanmean(CondB2.powspctrm(:,1:length(chans),:,tindx),1));
    assignin('base','meanCondA1',meanCondA1);
    assignin('base','meanCondB1',meanCondB1);
    assignin('base','meanCondA2',meanCondA2);
    assignin('base','meanCondB2',meanCondB2);
end

%% If we are interested in interaction effects (if size(DataIn,2) ==4) then find distance between 2 levels of the test


if size(cfg.DataIn,2) == 4
    fprintf('****Input cell array of length 4**** \n****Analysing Interaction Effects**** \n')
    meanCondA = meanCondA1 - meanCondA2;
    meanCondB = meanCondB1 - meanCondB2;
    
    mcfg = [];
    mcfg.parameter = 'powspctrm';
    mcfg.operation = 'subtract';
    
    CondA = ft_math(mcfg,CondA1, CondA2);
    CondB = ft_math(mcfg,CondB1, CondB2);
    assignin('base','CondA',CondA);
    assignin('base','CondB',CondB);
end

%% Calculate Non-Parametric Statistics on Time Frequency Data

if strcmp(Mode,'multielec')==1
    
    for fcnt=1:length(foi)   %Find the raw effect for each frequency of interest
        
        if length(size(meanCondA))>2
            lenf=length(CondA.freq);
            temp = ones(1,lenf).*foi(fcnt);
            [~,ind(fcnt)]=min(abs(CondA.freq-temp));
            raweffect{fcnt,1}=meanCondA(:,ind(fcnt),:)-meanCondB(:,ind(fcnt),:);
        else
            raweffect{fcnt,1} = meanCondA - meanCondB;
        end
        
    end
    
    assignin('base','raweffect',raweffect);
    
    layout_temp = '/Users/bolger/Documents/MATLAB/fieldtrip-master/template/neighbours/biosemi64_neighb.mat';
    cfg=[];
    cfg.method= 'triangulation';
    %cfg.template = layout_temp;
    %cfg.neighbourdist = 0.25;
    cfg.layout = layout;
    neighbs=ft_prepare_neighbours(cfg,CondA);
    assignin('base','neighbs',neighbs);
    
    
    len_foi = size(foi,1);
    statsout=cell(len_foi,1);  %define variables
    
    for count=1:len_foi
        count
        TFstatcfg=[];
        TFstatcfg.channel=chans;
        TFstatcfg.frequency =[foi(count,:)];
        TFstatcfg.latency=latency;
        TFstatcfg.avgovertime='no';
        TFstatcfg.avgoverchannel='no';
        TFstatcfg.avgoverfreq = freqavg;
        TFstatcfg.method = 'montecarlo';
        TFstatcfg.statistic = 'ft_statfun_depsamplesT';                          % activation vs baseline test
        TFstatcfg.correctm = 'cluster';                                               % method to deal with multiple comparisons
        TFstatcfg.clusteralpha = 0.05;
        TFstatcfg.clusterstatistic = 'maxsum';
        TFstatcfg.tail = 0;
        TFstatcfg.clustertail= 0;
        TFstatcfg.alpha=0.05;
        TFstatcfg.minnbchan = 2;            %the minimum number of neighbourhood channels required for selected channel to be included in the clustering algorithm
        TFstatcfg.numrandomization = 2000;  %number of draws from the permutation distribution
        TFstatcfg.neighbours = neighbs;
        
        Design=zeros(2,subnum*2);
        Design(1,1:subnum)=1:subnum;
        Design(1,subnum+1: subnum*2)=1:subnum;
        Design(2,1:subnum)=1;
        Design(2,subnum+1: subnum*2)=2;
        TFstatcfg.design=Design;
        TFstatcfg.uvar=1;
        TFstatcfg.ivar=2;
        
        statsout{count,1}=ft_freqstatistics(TFstatcfg,CondA,CondB);
        statsout{count,1}.raweffect = squeeze(raweffect{count,1});
        statsout{count,1}.logprob = -log10(statsout{count,1}.prob);
        assignin('base','statsout',statsout);
        
        hit=0;
        % ************Find significant positive clusters****************************
        if isfield(statsout{count,1},'posclusters')==1
            if isfield(statsout{count,1}.posclusters,'prob')==1
                PosSig=find([statsout{count,1}.posclusters.prob]<= TFstatcfg.alpha);
                if isempty(PosSig)==0
                    hit=hit+1;
                end
            else
                display('No Positive Clusters!')
                PosSig=[];
            end
        else
            PosSig=[];
        end
        % ***********Find significant negative clusters*****************************
        if isfield(statsout{count,1},'negclusters')==1
            if isfield(statsout{count,1}.negclusters,'prob')==1
                NegSig=find([statsout{count,1}.negclusters.prob]<= TFstatcfg.alpha);
                if isempty(NegSig)==0
                    hit=hit+1;
                end
            else
                display('No Negative Clusters!');
                NegSig=[];
            end
        else
            NegSig=[];
        end
        % If there are neither positive nor negative clusters.
        if isempty(PosSig)==1 && isempty(NegSig)==1
            display('No Significant Clusters so nothing to plot...sorry!!')
        end
        
        %% IF THERE ARE SIGNIFICANT CLUSTERS PLOT THE RESULTS ON TOPOGRAPHIC MAPS OVER TIME.
        
        if hit>0
            disp(horzcat('********************************Plotting for frequency ',num2str(foi(count)),'Hz***********************************'));
            
            assignin('base','time_steps',time_steps);
            
            rws=1; cols = length(time_steps);
            
            scrsz = get(groot,'ScreenSize');
            scrsz(1,4) = scrsz(1,4)/4;
            hndl=figure('Position',scrsz); set(hndl,'Color',[1 1 1]);
            
            for cnt = 1:size(time_steps,1)  % For each time window
                cnt
                
                
                %                 [axisl,axisb,axisw,axish] = CREx_prepare_subplot(hndl,cols,rws,cnt);   %Call of function to define subplot sizes
                %                 h=subplot('position',[axisl,axisb,axisw,axish]);
                
                subplot(rws, cols, cnt);
                
                StatOut2=[];
                raweffect_curr = [];
                
                Selcfg=[];     %Define the current latency
                Selcfg.latency = [time_steps(cnt,1) time_steps(cnt,2)];
                Selcfg.avgovertime='no';  %Do not average over the time window
                StatOut2 = ft_selectdata(Selcfg,statsout{count,1});
                assignin('base','StatOut2',StatOut2)
                
                toindx = [time_curr >=time_steps(cnt,1) & time_curr<=time_steps(cnt,2)];
                assignin('base','toindx',toindx);
                
                %Find if those electrodes that are significant over the
                %entirety of this time window
                sig_int = zeros(64,1);
                for chancnt = 1:64
                    sig_int(chancnt) = all(squeeze(StatOut2.mask(chancnt,:))); %Will return 1 if all time points in the current window are significant (all ones);
                end
                assignin('base','sig_int',sig_int);
                
                % Average the raw effect over current window of interest.
                raweffect_curr = nanmean(raweffect{count,1}(:,toindx),2);
                StatOut2.raweffect_curr = raweffect_curr;
                % Find the mask for the current time window
                mask_new = squeeze(statsout{count,1}.mask);
                mask_mean = sum(mask_new(:,toindx),2);
                mask_mean(mask_mean>0) = 1;
                StatOut2.mask = mask_mean;
                
                % Find the electrodes to highlight.
                A = cell(size(StatOut2.mask,2),1);
                for tcnt=1:size(StatOut2.mask,2)
                    Einds=arrayfun(@(x) find([x.prob(:,:,tcnt)<=min(StatOut2.prob)]),StatOut2,'UniformOutput',false);
                    A{tcnt,1} = Einds{1,1}';
                end
                
                
                prob_curr = StatOut2.prob;
                max_prob = min(prob_curr,[],2);
                StatOut2.logprob = -10*log10(max_prob);
                
     
                    topocfg=[];
                    topocfg.parameter= 'logprob';   %'';logprob-log10(.025)
                    topocfg.zlim = [0 3];
                    %topocfg.xlim=[time_steps(cnt,1) time_steps(cnt,2)];
                    topocfg.layout = layout;
                    topocfg.baseline='no';
                    topocfg.marker='off';
                    topocfg.markersymbol='.';
                    topocfg.markersize=24;
                    topocfg.comment='xlim';
                    topocfg.commentpos='title';
                    topocfg.highlight='on';
                    topocfg.highlightsymbol='.';
                    topocfg.highlightsize = 16;
                    
                    if sum(sig_int)>0
                        topocfg.highlightchannel = find(sig_int);
                    end
                    
                    topocfg.highlightcolor = [0 0 0];
                    topocfg.shading = 'interp';
                    topocfg.style = 'fill';
                    topocfg.colormap = colormap('parula');
                    
                    
                    ft_topoplotTFR(topocfg,StatOut2);
  
                
            end
            
        else
            disp(horzcat('********************************Nothing to plot for frequency ',num2str(foi(count)),'Hz***********************************'));
        end
        
    end
    
    %%
    
elseif strcmp(Mode,'indivelec')==1
    %% CARRY OUT CLUSTER-BASED PERMUTATION TEST ON INDIVIDUAL ELECTRODE DATA
    length(chans)
    statsout=cell(1,length(chans));
    raweffect=cell(1,length(chans));
    f1=figure; set(f1,'Color',[1 1 1],'Position',[100 100 1400 1000]);
    s=zeros(length(chans),1);
    chan_alpha = zeros(length(chans),1);
    
    for ecnt=1:length(chans)
        
        cfg=[];
        cfg.method='montecarlo';
        cfg.statistic= 'ft_statfun_depsamplesT';  %unit of observation here is subject or time point (if time point means that need to use indepsamplesT)
        cfg.correctm='cluster';
        cfg.clusteralpha=0.05;
        cfg.clustertail= 0;
        cfg.clusterstatistic='maxsum';
        cfg.tail= 0; % two-sided test
        cfg.correcttail='prob';
        cfg.alpha=0.05;
        cfg.numrandomization = 2000;
        cfg.latency=latency;
        cfg.frequency = [foi(1) foi(end)];
        
        design=zeros(2,subnum*2);  %initialise the design matrix
        design(1,1:subnum)=1:subnum;
        design(1,subnum+1:subnum*2)=1:subnum;
        design(2,1:subnum)=1;
        design(2,subnum+1:subnum*2)=2;
        
        cfg.design=design;
        cfg.ivar=2;
        cfg.uvar=1;
        cfg.neighbours=[];
        cfg.channel=chans{1,ecnt};   %define the current electrode of interest
        
        chanind=find(strcmp(chans{1,ecnt},CondA.label));
        statsout{1,ecnt}=ft_freqstatistics(cfg,CondA,CondB);
        
        time=CondA.time;
        ti=find(time>=cfg.latency(1) & time<=cfg.latency(2));
        assignin('base','ti',ti);
        ti = cat(2,ti(1)-1,ti);
        
        
        meanCondA_pow=squeeze(mean(CondA.powspctrm(:,:,:,ti),1)); disp(size(meanCondA_pow));
        meanCondB_pow=squeeze(mean(CondB.powspctrm(:,:,:,ti),1)); disp(size(meanCondB_pow));
        assignin('base','meanCondB_pow',meanCondB_pow);
        assignin('base','meanCondA_pow',meanCondA_pow);
        raweffect{1,ecnt}=meanCondA_pow(chanind,:,:)-meanCondB_pow(chanind,:,:);
        S=raweffect{1,ecnt}; %freq X time
        statsout{1,ecnt}.effect=S(:,:,1:end); %may need to subtract 1
        statsout{1,ecnt}.logprob = -log10(statsout{1,ecnt}.prob);
        
        [x,y,z]=find(isnan(statsout{1,ecnt}.effect));
        
        statsout{1,ecnt}.effect(x,y,z)=0;
        assignin('base','statsout_curr', statsout{1,ecnt});
        
        %PLOT TIME-FREQUENCY PLOT WITH SIGNIFICANT TIME-FREQUENCY REGIONS
        %HIGHLIGHTED
        figcfg=[];
        figcfg.channel = chans{1,ecnt};
        figcfg.baseline = 'no';
        figcfg.renderer = 'openGL';
        figcfg.color = 'yes';
        figcfg.parameter = 'logprob';
        figcfg.maskstyle = 'opacity';
        figcfg.maskparameter = 'mask';
        figcfg.masknans = 'yes';
        figcfg.maskalpha = 0.6;
        figcfg.zlim= [0 3];
        figcfg.xlim = latency;
        figcfg.layout=layout;
        
        rws=3; cols=ceil(length(chans)/rws);
        subplot(rws,cols,ecnt);
        ft_singleplotTFR(figcfg,statsout{1,ecnt});
        title(horzcat('Electrode: ',figcfg.channel));
        s(ecnt)=gca;
        set(s(ecnt),'HitTest','on','SelectionHighlight','on','UserData',{figcfg,statsout{1,ecnt}},'ButtonDownFcn',@plotsingleTF_elec);    %call of function
        
        chan_alpha(ecnt)= min(unique(min(statsout{1,ecnt}.prob)));
    end
    
    assignin('base','statsout', statsout);
    assignin('base','chan_alpha',chan_alpha);
    
end

end
