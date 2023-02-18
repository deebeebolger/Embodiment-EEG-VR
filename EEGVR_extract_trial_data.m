
function filesave = EEGVR_extract_trial_data(subjcurr,dirbase)
% Date: March 2018                  Programmed by: D. Bolger
% Function to bring together all the trial-level files output by Unity (VR
% system) to create a single *.txt file containing trial-level information
% concerning the following:
% 1. Trial index.
% 2. Auditory stimulus (verb) presented on this trial.
% 3. ObjTarget (not used)
% 4. ObjDistractor (recorded by not used).
% 5. Trialtype - "go" or "no-go".
% 6. Presentation side - right or left.
% 7. Bad trial - specifying if subject-response on this trial was correct or incorrect.
% Output : filesave : the title of the stimulus data file for the current
% subject.
%-------------------------------------------------------------------------------------------

dircurr = strcat(dirbase,'StimData/Trials/'); % Path to Unity trial-level files. 
X = dir(dircurr); 
filenoms = {X(~[X.isdir]).name};
filenoms =  filenoms(~strcmp(filenoms,'.DS_Store'));
trialnum = length(filenoms);
stim_save = strcat(dirbase,'StimData/');

% Prepare stimulus data *.txt file. 
filesave = strcat(subjcurr,'_stimdata.txt'); 
fid1 = fopen(strcat(stim_save,filesave),'a+'); 
fprintf(fid1,'//TrialNum Verb ObjTarget ObjDistractor TrialType PresSide BadTrial');

% Run for each trial. 
for tcnt = 1:trialnum
    
    display(tcnt)
    
    tagcurr = strcat('Trial_Order',num2str(tcnt-1),'_'); 
    I = strfind(filenoms,tagcurr); 
    indx = find(cellfun(@isempty,I)==0); 
    
    % Open the file at index indx and read in first 5 lines
    fid3 = fopen(strcat(dircurr,filenoms{1,indx}));
    condcurr = textscan(fid3, '%s%s','Delimiter','\t');   % This should be a 1 x 5 cell array.
    fclose(fid3);
    clear fid3; 
    
    verbcurr = [];
    objtarget = [];
    objdist = [];
    go_nogo =[];
    objside = []; 
    
    verbcurr = condcurr{1,2}{1};
    objtarget = condcurr{1,2}{2};
    objdist = condcurr{1,2}{3};
    go_nogo = condcurr{1,2}{4};
    objside = condcurr{1,2}{5}; 
    
    fprintf(fid1,'\n%d %s %s %s %s %s',tcnt, verbcurr,objtarget,objdist,go_nogo,objside);
    
    if ~isempty(strfind(filenoms{1,indx},'REPEAT_AT_THE_END_OF_BLOCK'))
        badtrial = 'badtrial';
        fprintf(fid1,' %s',badtrial); 
    else
        badtrial = 'goodtrial';
        fprintf(fid1,' %s',badtrial);
    end
    
    
end

fclose(fid1); 



