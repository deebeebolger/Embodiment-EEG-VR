function [tdata_extend, tdata_lims,T]=CREx_trialextend_mirror(tdata,catopt,fs,time)
%Date: 16-05-2017                Programmed by: D. Bolger
%function to concatenate data onto both the beginning and end of a single
%trial of data by creating a mirror image of the data. 
%Input data: 
% - tdata : single trial data for a single electrode
% - catopt: 'start' 'end' or 'both'
% - fs: sampling rate
% - time: time vector of input data
% Output:
% - tdata_extend: extended version of input data
% - tdata_lims: indices of start and end of the original input data
%***************************************************************
len=size(tdata,2); 
tnum=size(tdata,3);
indt0=find(time==0);    %Find the T0 indice

if isempty(indt0)    %To deal with lack of 0 point after downsampling
    idx=nearest(time,0);
    time=time-time(idx);
    indt0=find(time==0); 
end 
repnum=[3 3]; 

tdata_flipend=flip(tdata(:,indt0:end-1,:)); 
addtotalend=repmat(tdata_flipend,[1,repnum(1)]);                          

tdata_flipstart=flip(tdata(:,1:indt0,:));
addtotalstart=repmat(tdata_flipstart,[1,repnum(2)]);

assignin('base','addtotalstart',addtotalstart);
assignin('base','addtotalend',addtotalend); 
assignin('base','tdata',tdata);
assignin('base','indt0',indt0); 

if strcmp(catopt,'start')
    tdata_extend=cat(2,addtotalstart,tdata);
    tdata_lims=[length(addtotalstart)+1 (length(addtotalstart)+1)+len];
    Tzero=indt0+length(addtotalstart); 
    T= cat(2,flip((1/fs.*(1:Tzero-1)).*-1),time(indt0:end)); 
elseif strcmp(catopt,'end')
    tdata_extend=cat(2,tdata,addtotalend);
    tdata_lims=[1 len];
    len2=len-indt0; 
    tadd=1/fs.*(1:(len2+size(addtotalend,2)));
    T= cat(2,time(1:indt0),tadd); 
elseif strcmp(catopt,'both')
    lowend=cat(2,addtotalstart,tdata);
    tdata_extend=cat(2,lowend,addtotalend);
    tdata_lims=[size(addtotalstart,2)+1 (size(addtotalstart,2)+1)+len]; 
    Tzero=indt0+length(addtotalstart); 
    len2=len-indt0; 
    tadd=1/fs.*(1:len2+size(addtotalend,2));
    T= cat(2,flip((1/fs.*(1:Tzero-1)).*-1),0,tadd);
end

 

end % end of main function