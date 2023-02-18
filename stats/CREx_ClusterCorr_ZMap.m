function CREx_ClusterCorr_ZMap(tfinput,chanx,axhdl,BLine,fhndl)
% Function to test the change in time-frequency power for a single
% condition compared to a baseline interval. It applies nonparametric
% permutation testing and corrects for multiple comparisons at the
% cluster-level. 
% The synchronisation/desynchronisation is expressed as z-score.
% The plot, therefore, is a zmap.
% This works for a single electrode.
% The idea behind this test is to test if the time-freuency power of itself
% is statistically significant. 
% Input:
% tfinput: Structure with fields - powspctrm, freq, time, labels
%          structure of powspctrm array: trials/subjects X chan X freq X time
% chanx: current channel of interest
%**********************************************************************
F = tfinput.freq;  %the frequency vector
foi = [4 35];
findx1 = [dsearchn(F',foi(1)):dsearchn(F',foi(2))];
assignin('base','findx1',findx1);
eegpower=tfinput.powspctrm(:,:,findx1,:);                     %the calculated frequency spectrum
freqx = F(findx1);                        
num_freqx=size(eegpower,3);                     %the number of frequency components 
T=tfinput.time;                                 %the time vector
indz = dsearchn(T',0);
tlims = [T(indz) T(end)];                                %the time limits of the epoch
tindx=[T>=tlims(1) & T<=tlims(2)];              %Find the time points falling within epoch limits 
tftimes_all = T(tindx); 
ntrials = size(eegpower,1);                       %Number of trials/subjects
chanindx=find(strcmp(tfinput.label,chanx));     %Only want the 64 electrodes
eegpow=squeeze(eegpower(:,chanindx,:,tindx));   %Single Electrode: the power spectrum array has format: trials X freq X time
assignin('base','eegpow',eegpow);


bl_lims = BLine;   %Specify the baseline limits
tpoints = numel(tftimes_all);                     %Number of time samples
voxel_pval   = 0.1;                            %Uncorrected map p-value
cluster_pval = 0.1;                            %Cluster-correction p-value

%Number of permutations
nperms = 1000;

blidx(1) = dsearchn(tftimes_all',bl_lims(1)); %Find the nearest time-point to the baseline time-limits
blidx(2) = dsearchn(tftimes_all',bl_lims(2));
assignin('base','blidx',blidx)
assignin('base','bl_lims',bl_lims)


% compute t-test of baseline-poststimulus difference
realbls = squeeze(mean(eegpow(:,:,blidx(1):blidx(2)),3));                          %Find the mean power over the baseline.
realavg = squeeze(10*log10(bsxfun(@rdivide, mean(eegpow,1), mean(realbls,1))));    %Calculate the baseline normalised log power
assignin('base','realbls',realbls);
assignin('base','realavg',realavg);

% initialize matrices (null hypothesis) 
permuted_maxvals = zeros(nperms,2,num_freqx);
valperm    = zeros(nperms,num_freqx,tpoints);
max_clust_info   = zeros(nperms,1);

for pcnt=1:nperms
    cutpoint = randsample(2:tpoints-diff(blidx)-2,1);    %Randomly select a single cutoff point.
    valperm(pcnt,:,:) = 10*log10(bsxfun(@rdivide,mean(eegpow(:,:,[cutpoint:end 1:cutpoint-1]),1),mean(realbls,1))); %computing the decibel change from baseline at each iteration of permutation testing
end

zmap = (realavg-squeeze(mean(valperm))) ./ squeeze(std(valperm));  %Express as z-score
threshmean = realavg;
threshmean(zmap < norminv(1-voxel_pval)) = 0;

assignin('base','zmap',zmap);


%% CARRY OUT CLUSTER CORRECTION ON PERMUTED DATA
maxclust_info=zeros(length(nperms),1); 

for permi = 1:nperms
    
    % for cluster correction, apply uncorrected threshold first to get maximum cluster sizes
    fakecorrsz = squeeze((valperm(permi,:,:)-mean(valperm,1)) ./ std(valperm,[],1) );
    assignin('base','fakecorrsz',fakecorrsz);
    assignin('base','valperm',valperm)
    fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
    
    % get number of elements in largest above-threshold cluster
    cluster_info = bwconncomp(fakecorrsz);                                       % scan the thresholded t-values for clusters
    maxclust_info(permi) = max([ 0 cellfun(@numel,cluster_info.PixelIdxList) ]); % find the maximum cluster sizes
   
end

assignin('base','cluster_info',cluster_info);
assignin('base','maxclust_info',maxclust_info);

zmap_corr = zmap;                                               % apply cluster-level corrected threshold
zmap_corr(zmap_corr<norminv(1-voxel_pval))=0;            % uncorrected pixel-level threshold
cluster_info = bwconncomp(zmap_corr);                           % find isolated clusters and remove those smaller than cluster size threshold
clusterinfo = cellfun(@numel,cluster_info.PixelIdxList);       %
cluster_thresh = prctile(maxclust_info,100-cluster_pval*100);   % 95% percentile as cluster_pval = 0.05
cluster_rm = find(clusterinfo<cluster_thresh);                 % clusters falling out the threshold are removed

% remove clusters
for i=1:length(cluster_rm)
    zmap_corr(cluster_info.PixelIdxList{cluster_rm(i)})=0;
end

assignin('base','zmap_corr',zmap_corr);
assignin('base','tlims',tlims); 
assignin('base','zmap',zmap);

%% PLOT THE CLUSTER CORRECTED Z-MAP RESULTS

contourf(tftimes_all.*1000,freqx,zmap,40,'linecolor','none')
set(axhdl,'clim',[-3 3],'xlim',tlims.*1000)
title(chanx)

set(fhndl,'CurrentAxes',axhdl);
set(axhdl,'HitTest','on','SelectionHighlight','on','UserData',{tftimes_all,freqx,zmap_corr,chanx,tlims});
set(axhdl,'ButtonDownFcn',@plotsingle_map)


end % function end 

function plotsingle_map(hdl,~)

 D=get(hdl,'UserData');
 T=D{1,1};
 freqs=D{1,2};
 zm_corr=D{1,3};
 chanoi=D{1,4}; 
 tlim=D{1,5}; 
 
f1=figure; set(f1,'PaperUnits','normalized','PaperPosition',[0.0235308 0.0272775 0.894169 0.909249],'Color',[1 1 1]);
orient portrait; axis ('normal');
 
 contourf(T.*1000,freqs,zm_corr,40,'linecolor','none');
 set(gca,'clim',[-3 3],'xlim',tlim.*1000);
 title(strcat('Cluster-corrected Zmap: ',char(chanoi)));
 
end 
