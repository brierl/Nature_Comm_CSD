%% WRAPPER: HGB, GCAMP, EEG DATA PROCESSING for CSDs

% This wrapper runs through the image and EEG processing steps for a mouse.
% This first section contains information for locating and loading the
% correct mouse's data. rawFolder is the data main directory. Raw data is 
% sorted into rawFolder/Date/ directories. Processed data is sorted into 
% rawFolder/Analysis/Date/ directories. Assumes naming convention 
% Date-Mouse-trials{i} for mouse image stack is used, 
% i loops through runs, below. 

% Wrapper calls 3 sub-routines, makemask.m,
% proc_Imaging.m and proc_EEG.m. 

clear
close all

%data paths
rawFolder='/Users/lindsey/Documents/COVID-19/CSD/'; %path to main raw data directory
foldOut=[rawFolder,'Analysis/']; %directory to output processed data

%mouse info
Date='170915'; %date mouse was imaged
gcamp_file=1; %flag for if mouse has gcamp, 1=yes, 0=no
EEG_file=1; %flag for if mouse has eeg, 1=yes, 0=no
foldDate=[rawFolder,Date,'/']; %organized into folders by date, for raw data
files = dir([foldDate '*.tif']); %extract names of all files
files = {files.name}';
if exist([foldOut,Date], 'file') == 0
    mkdir([foldOut,Date]) %organized into folders by date, for saved processed data
end
Mouse='GCAMPMs1'; %mouse name
trials={'kcl1.tif','kcl2.tif'}; %run names

%% Analysis begins
    
% Make masks
if exist([foldOut,Date,'/',Date,'-',char(Mouse),'-mask.tif'])>0
else
    [mask]=makemask(Mouse,Date,foldOut,trials,rawFolder);
end
% mask: .tif file with same x,y dimensions as raw data setting all
% non-brain regions to 0. Saves mask as Date-Mouse-mask.tif

for trial=1:length(trials) %loop through runs    
              
    filename=[Date,'-',char(Mouse),'-',char(trials{trial})]; %assemble filename using naming convention     
    disp(['T R I A L :   ',filename,])
    if exist([foldOut,Date,'/',Date,'-',char(Mouse),'-',char(trials{trial}(1:end-4)),'-DATA.mat'])>0
    else
        if trial==1 %initiate to empty if first trial of mouse, will process 1st run now...
            rs_gcamp6corr=[];
            rs_colorchannels=[];
            rs_gcamp6u=[];
        else
            %load baseline resting state data aka first run of mouse. Will
            %use this baseline data to process CSD runs now...
            load([foldOut,Date,'/',Date,'-',char(Mouse),'-',char(trials{1}(1:end-4)),'-DATA.mat'],'color_channels','gcamp6corr','gcamp6u');
            rs_gcamp6corr=gcamp6corr; %ratio corrected gcamp for hgb confound (Wright et al., 2017). pixel x pixel x frame
            rs_colorchannels=color_channels; %individual LED channels for oximetry. pixel x pixel x LED# x frames
            rs_gcamp6u=gcamp6u; %uncorrected gcamp for hgb confound. pixel x pixel x frame
        end
        
        %Process imaging data
        [data_dot,gcamp6corr,WL,info,color_channels,gcamp6u]=proc_Imaging(rawFolder,filename,Date,gcamp_file,rs_gcamp6corr,rs_gcamp6u,rs_colorchannels);
        % data_dot=oxy(:,:,1,:) or deoxy(:,:,2,:) hgb. pixel x pixel x contrast x frame
        % gcamp6corr=ratio corrected gcamp for hgb confound (Wright et al., 2017). pixel x pixel x frame
        % WL=white light image. pixel x pixel
        % color_channels=individual LED channels for oximetry. pixel x pixel x LED# x frames
        % gcamp6u=uncorrected gcamp for hgb confound. pixel x pixel x frame
        save([foldOut,Date,'/',Date,'-',char(Mouse),'-',char(trials{trial}(1:end-4)),'-DATA.mat'],'-v7.3','data_dot','gcamp6corr','gcamp6u','WL','info','color_channels');
    
    end

end

if EEG_file==1 %if mouse had concurrent EEG
    [gcamp_eeg,fs]=proc_EEG(Mouse,Date,rawFolder,foldOut);
    %gcamp_eeg: run number x length of processed EEG trace. EEG trace with
    %   mean subtracted.
    %fs: sample rate for EEG, default 10000    close all
end
