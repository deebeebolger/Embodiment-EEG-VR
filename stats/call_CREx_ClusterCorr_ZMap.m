% Date: 05-2019                 Programmed: D. Bolger
% Figure to call the function to carry out baseline vs. post-stimulus
% comparision via a non-parametric test and a cluster-based correction.
% This is carried out on an electrode-by-electrode basis and after each
% calculation, the zmap is plotted. 

InputData = GATFdata_all{2,1};                                                                                                                                                                       
lenf = length(InputData.freq);
Chanoi={InputData.label{1:64}};
Time=InputData.time;
baseline = [-0.2 0];
InputData.freq = InputData.freq(1:lenf);
InputData.powspctrm = InputData.powspctrm(:,:,1:lenf,:);

%% PREPARE THE FIGURE TO PLOT THE ERPs WITH ELECTRODE CONFIGURATION

hndl = figure; set(hndl,'PaperUnits','normalized','PaperPosition',[0.0235308 0.0272775 0.894169 0.909249],'Color',[1 1 1]);
orient portrait; axis ('normal')

xvals=zeros(length(chaninfo),1);
yvals=zeros(length(chaninfo),1);
pwidth    = 0.75;     % 0.75, width and height of plot array on figure
pheight   = 0.75;
axwidth=0.04;
axheight=0.08;

%Identify all the non-empty channel
notmtchans=cellfun('isempty',{chaninfo.theta});
notmtchans=find(~notmtchans);


%Read in channel locations file
[elocs,titres,theta,rads,inds]=readlocs(chaninfo(notmtchans));
channoms = strvcat(chaninfo.labels);
Th=pi/180*theta;           %convert degrees to radians

%Convert from polar to cartesian
[ycart,xcart]=pol2cart(Th,rads);
xvals(notmtchans)=ycart;
yvals(notmtchans)=xcart;

%Find the positions of all the channels
allchans=length(chaninfo);
mtchans=setdiff(1:allchans,notmtchans);        %find the channels indices common to both
allchans_sqrt=floor(sqrt(allchans))+1;

for i=1:length(mtchans)
    
    xvals(mtchans(i))=0.7+0.2*floor((i-1)/allchans);
    yvals(mtchans(i))=-0.4+mod(i-1,allchans)/allchans;
    
end

channoms2=channoms(1:64,:);
xvals=xvals(1:64);
yvals=yvals(1:64);

if length(xvals) > 1
    if length(unique(xvals)) > 1
        gcapos = get(gca,'Position'); axis off;
        xvals = (xvals-mean([max(xvals) min(xvals)]))/(max(xvals)-min(xvals)); % this recenters
        xvals = gcapos(1)+gcapos(3)/2+pwidth*xvals;                                 %this controls width of plot
    end
end

gcapos = get(gca,'Position'); axis off;
yvals = gcapos(2)+gcapos(4)/2+pheight*yvals;  % controls height of plot
ho=zeros(length(Chanoi),1);
sig_elecs=zeros(length(Chanoi),length(Time));
sig_times=cell(length(Chanoi),1);
nosigs=0;

Axes = cell(length(Chanoi),1);

for cnt=1:length(Chanoi)
    
    disp(Chanoi{1,cnt}) 
    xcenter=xvals(cnt);
    ycenter=yvals(cnt);
    Axes{cnt,1} = [Axes{cnt,1} axes('Units','normalized','Position', [ycenter-axheight/2 xcenter-axwidth/2 axheight axwidth])];
    
    CREx_ClusterCorr_ZMap(InputData,Chanoi{1,cnt},Axes{cnt,1},baseline,hndl);
    hold on;
end

