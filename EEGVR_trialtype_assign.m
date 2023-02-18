function [EEG,condcurr] = EEGVR_trialtype_assign(EEG,sujnum,dirbase)
%% Function to assign the trial-type, "go", "nogo" and "badtrial" to the events.
% Date: 01-05-2018         Programmed by: D. Bolger
% Input Data:
% EEG - EEG structure of the current subject.
% sujnum - the title (number and list letter) of the current subject.
% dirbase - the current base directory for the current subject.
% Output Data:
% EEG - EEG structure of the current subject with events trial-types added
% to event field.
% condcurr - cell array containing data from the StimData file.
% Use as: [EEG,condcurr] = EEGVR_trialtype_assign(EEG,sujnum,dirbase)
%******************************************************************************

stimfile = strcat(sujnum,'_stimdata.txt');
stimdir = fullfile(dirbase,'StimData',stimfile);

%% OPEN STIMULUS TEXT FILE FOR CURRENT SUBJECT.

fid = [];
fid = fopen(stimdir);
condcurr = textscan(fid, '%d%s%s%s%s%s%s','CommentStyle','//');   % This should be a 1 x 6 cell array.
fclose(fid);

verbs = unique(condcurr{1,2});
goodtrials = condcurr{1,7}; 
events_all = {EEG.event.type};

for counter = 1:length(verbs)
    
    % Find the indices of all instances of the current verb in the stimulus presentation
    % sequence (condcurr)
    vindx = strcmp(verbs{counter,1},condcurr{1,2});  % The length of vindx and x should be the same. 
    gdtrial = goodtrials(vindx);
    
    x = find(strcmp(verbs{counter,1},events_all));
    
    for counter2 = 1:length(x)
        if strcmp(events_all{1,x(counter2)-1},'go')  && strcmp(gdtrial{counter2,1},'goodtrial')
            events_all{1,x(counter2)} = strcat(events_all{1,x(counter2)},'-go');
            EEG.event(x(counter2)).type = events_all{1,x(counter2)} ;
        elseif strcmp(events_all{1,x(counter2)-1},'nogo')  && strcmp(gdtrial{counter2,1},'goodtrial')
            events_all{1,x(counter2)} = strcat(events_all{1,x(counter2)},'-nogo');
            EEG.event(x(counter2)).type = events_all{1,x(counter2)} ;
        elseif strcmp(events_all{1,x(counter2)-1},'go') && strcmp(gdtrial{counter2,1},'badtrial')
            events_all{1,x(counter2)} = strcat(events_all{1,x(counter2)},'-bad');
            EEG.event(x(counter2)).type = events_all{1,x(counter2)};
        elseif strcmp(events_all{1,x(counter2)-1},'nogo') && strcmp(gdtrial{counter2,1},'badtrial')
            events_all{1,x(counter2)} = strcat(events_all{1,x(counter2)},'-bad');
            EEG.event(x(counter2)).type = events_all{1,x(counter2)};
        end 
    end
end


end %end of function 

