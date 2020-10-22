function [gcamp_eeg,fs]=proc_EEG(Mouse,Date,rawFolder,foldOut)
% proc_EEG.m loads an EEG array for a particular mouse, subtracts
% the mean EEG trace for each run (trial), plots the result and saves
% the plot and processed EEG. Assumes naming convention 
% Date-Mouse-EEG.mat for loading EEG for a mouse. Saves processed EEG as
% Date-Mouse-EEGprocessed.mat and plot as Date-Mouse-EEGprocessed.jpg

% IN: 
%   Mouse: mouse name
%   Date: date mouse was imaged
%   rawFolder: rawFolder/Date/ is directory containing raw data
%   foldOut: foldOut/Date/ is directory to save processed EEG and plot to

% OUT:
%   gcamp_eeg: run number x length of processed EEG trace. EEG trace with
%       mean subtracted.
%   fs: sample rate for EEG, default 10000

    run=[Date,'-',Mouse,'-EEG.mat']; %raw EEG file name
    runout=[Date,'-',Mouse,'-EEGprocessed.mat']; %processed EEG save name
    
    eegdata=cell2mat(struct2cell(load([rawFolder,Date,'/',run],'data'))); %load EEG matrix
    eegdata=[eegdata zeros(1,60000000-length(eegdata))]; %make even, assuming each eeg run is 10 minutes, fs 10000
    eegstart=cell2mat(struct2cell(load([rawFolder,Date,'/',run],'datastart'))); %load flag for beginning of recording
    eegend=cell2mat(struct2cell(load([rawFolder,Date,'/',run],'dataend'))); %load flag for end of recording
    fs=cell2mat(struct2cell(load([rawFolder,Date,'/',run],'samplerate'))); %load sample rate, default 10000
    fs=fs(1);

    trial_length=min(eegend-eegstart); %find the lower bound of length that constitutes a trial 
    gcamp_eeg=zeros(length(eegstart),6000000); %assumes each eeg trial is 10 minutes, fs 10000, initialize matrix for processed EEG
    for n=1:length(eegstart) %proc per run
    	
        IndividualRun=eegdata(eegstart(1,n):eegend(1,n)); %isolate individual run
        if length(IndividualRun)<trial_length %false start, not a full run, try next flag that signifies beginning of recording
        	continue
        else
            runs=IndividualRun(5:trial_length)-mean(IndividualRun(5:trial_length));  %subtract the mean, trim off first few samples and end so all trials have same length
        end
        gcamp_eeg(n,1:length(runs))=runs; %load into matrix for processed EEG
        
    end
    
    save([foldOut,Date,'/',runout], 'gcamp_eeg','fs');
    %% Plot result
       
	figure('Position',[0 0 1500 300]);
    hold on;
    for i = 1:size(gcamp_eeg,1) %plot each processed EEG run, yellow vertical lines indicate trial divisions
        subplot('Position',[(i-1)/size(gcamp_eeg,1) 0 1/size(gcamp_eeg,1) 1]);
        hold on;
        plot(gcamp_eeg(i,:)) ;
        plot([0 0],[-0.005 0.005],'-y');
        ax = gca;
        ax.Visible = 'off';
        ylim([-0.005 0.005]);
        xlim([0 trial_length]);
        title(i);
    end
    
    figTag=[foldOut,Date,'/',Date,'-',Mouse,'-EEGprocessed.jpg']; %save figure
    print ('-djpeg', '-r300', figTag);

end