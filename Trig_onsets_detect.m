function trialtrig =Trig_onsets_detect(DIn,time,fs,dircurr)
% Programmed by: Deirdre Bolger         Date: 10-01-2018
% This programme detects the onsets of the photodiode triggers. 
% It also detects the offsets. 
% Input data:
% DIn = channel data comprising photodiode trigger data (EEG.data(73,:))
% time = time vector (EEG.times)
% dircurr = the current subject-dependent path.
% sujlabel = the current subject label.
%*********************************************************************
%% LOAD IN THE PARAMETERS FILE DATA FOR CONDITION DEPENDENT PHOTODIODE DURATIONS
% outParameters.txt

fdir = strcat(dircurr,'StimData/outParameters.txt');
fid = fopen(fdir);
params = textscan(fid, '%s%s','Delimiter','\t');   % This should be a 1 x 5 cell array.
fclose(fid);

igo = ismember(params{1,1},'TimeRectGo');       % find duration for the go trials
godur = str2double(cell2mat(params{1,2}(igo)))*1000; 
inogo = ismember(params{1,1},'TimeRectNoGo');   % find duration for the no-go trials
nogodur = str2double(cell2mat(params{1,2}(inogo)))*1000;
ileft = ismember(params{1,1},'TimeRectLeft');   % find duration of left-presentation trig
leftdur = str2double(cell2mat(params{1,2}(ileft)))*1000; 
iright = ismember(params{1,1},'TimeRectRight'); % find duration of right-presentation trig.
rightdur = str2double(cell2mat(params{1,2}(iright)))*1000; 

%% **************************Extract the photodiode onsets****************************
D_raw = DIn(1,:);
D = detrend(D_raw,0);
thresh_val = 0;
D = D-thresh_val;

[b,a] = butter(2,8./(fs/2));
Dfilt = filtfilt(b,a,double(D));

%Half wave rectify the filtered signal and invert.
D_hwr=zeros(size(Dfilt));
ipos=find(Dfilt<0);
D_hwr(ipos)=Dfilt(ipos).*-1;

% Set all activity <= mean activity to zero.
D_hwr(D_hwr<=mean(D_hwr))=0; 

% Turn all trigger signals into step functions. 
D_hwr(D_hwr>0) = 1;

% Find onsets (diff(D_hwr == 1)) and offsets (diff(D_hwr == -1))
Dhwr_diff = diff(D_hwr);  % 1OD of D_hwr
Dhwr_diff = cat(2,Dhwr_diff,0);

[pks,minima,locs_pks,locs_min]= CREx_peakfinder(Dhwr_diff); 
assignin('base','locs_pks',locs_pks);
assignin('base','locs_min',locs_min);
assignin('base','Dhwr_diff',Dhwr_diff);

% Maybe include error message if locs_min and locs_pks do not have same
% length. 

onoffsets = [time(locs_pks);time(locs_min)]'; %locs_min(2:end)
durs = onoffsets(:,2) - onoffsets(:,1); 
onoffsets = cat(2,onoffsets,durs); 

onsets_all = nan(size(time));
offsets_all = nan(size(time));
onsets_all(Dhwr_diff==pks(1)) = 0;
offsets_all(Dhwr_diff==minima(1)) = 0;

figure;
subplot(2,1,1)
plot(time,D_hwr);
hold on
plot(time,onsets_all,'or','MarkerFaceColor','r')
hold on
plot(time,offsets_all,'og','MarkerFaceColor','g'); 
set(gca,'YLim',[0 1.5])
title('Photodiode Signal as Step-function: onsets and offsets'); 
subplot(2,1,2)
plot(time,D)
hold on
plot(time,onsets_all,'or','MarkerFaceColor','r');
hold on
plot(time,offsets_all,'og','MarkerFaceColor','g');
title('Original photodiode signal with onsets (red) and offsets (green)');

%% ****************Categorize the photodiode signal************************************

durlims = [nogodur godur leftdur rightdur;nogodur+150 godur+150 leftdur+500 rightdur+400]';
ttypes = {'nogo' 'go' 'left' 'right'}';
trialtrig = cell(length(durs),4); 

for Counter = 1 : 4
    indx = find(durs<=durlims(Counter,2) & durs>=durlims(Counter,1));
    for ncnt = 1:length(indx) 
        trialtrig{indx(ncnt),1} = ttypes{Counter};
        trialtrig{indx(ncnt),2} = onoffsets(indx(ncnt),1)./1000; 
        trialtrig{indx(ncnt),3} = durs(indx(ncnt))./1000; 
        trialtrig{indx(ncnt),4} = locs_pks(indx(ncnt)); 
    end
end



end 





