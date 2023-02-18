function [EEG] = EEGVR_addcovariates(EEG,verbcurr,Dirxls)
%% Date: 24-5-2018                Programmed by: D. Bolger
% Function to extract covariates such as word frequency and uniqueness
% point data from an excel file to integrate it into the subject-level
% data; the EEG structure. 
% The input data is the epoched data. 
% Input data: EEG structure of current subject with epoched data.
%             verbcurr - cell array of verbs found in current dataset.
%             Dirxls - path to xls file with word frequency data. 
% Output data: EEG structure with the covariates integrated into the "event" field of the structure.  
%-----------------------------------------------------------------

% Read line 1 for covariate type.
% Read lines 1 to 17 and columns A to AH

[~,~,alldata] = xlsread(Dirxls,'A1:AH17');
types = alldata(1,:);

wfreqlemfilm = alldata(:,ismember(types,'freqlemfilms2'));
wfreqlemfilm = wfreqlemfilm(2:end,1);
wfreqfilm = alldata(:,ismember(types,'freqfilms2'));
wfreqfilm = wfreqfilm(2:end,1);
verbs = alldata(:,ismember(types,'ortho'));
verbs = verbs(2:end,1);
puphon = alldata(:,ismember(types,'puphon'));
puphon = puphon(2:end,1);

%% For each verb, find instances of this verb in the current event field.
% If the verb exists, then extracts its freqlemfilm, freqfilm and puphon
% data and add them to the event field of the EEG structure. 

freqlemfilm = zeros(length(verbcurr),1);
freqfilm = zeros(length(verbcurr),1);
pointUphon = zeros(length(verbcurr),1);

for i = 1:length(verbcurr)

    indx = strncmp(verbs,verbcurr{1,i},3);

    freqlemfilm(i) = str2double(wfreqlemfilm{find(indx),1});
    freqfilm(i) = str2double(wfreqfilm{find(indx),1});
    pointUphon(i) = puphon{find(indx),1};

    EEG.event(i).UPphon = pointUphon(i,1);
    EEG.event(i).FreqLemFilm = freqlemfilm(i,1);
    EEG.event(i).FreqFilm = freqfilm(i,1);
    

end


end % end of function main. 

