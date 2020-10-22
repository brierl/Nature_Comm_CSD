function [mask]=makemask(Mouse,Date,foldOut,trials,rawFolder)
% makemask.m creates a .tif file the same dimensions as each .tif file in
% a mouse image stack. Uses roipoly.m so user can click around brain
% regions. Once the brain has been traced, right click--> create mask will
% set all non-brain regions to 0. Assumes naming convention 
% Date-Mouse-trials{i} for mouse imaging is used, where
% i loops through runs. Saves mask as Date-Mouse-mask.tif

% IN: 
%   Mouse: mouse name
%   Date: date mouse was imaged
%   foldOut: foldOut/Date/ is directory to save mask to
%   trials: run names-makes mask based on first trial for mouse
%   rawFolder: rawFolder/Date/ is directory containing raw data

% OUT:
%   mask: .tif file with same x,y dimensions as raw data setting all
%   non-brain regions to 0.

    disp(['Making mask for ',foldOut,Date,'/',Date,'-',char(Mouse)]);
    
    fcrun=[rawFolder,Date,'/',Date,'-',char(Mouse),'-',char(trials{1}),];
    greenframe=double(imread(fcrun,6)); %load green frame on first run for mouse
    greenframeall=greenframe./max(greenframe(:)); %normalize frame

    mask=double(roipoly(greenframeall)); %Click around brain regions, right click--> create mask
    imwrite(mask,[foldOut,Date,'/',Date,'-',char(Mouse),'-mask.tif'],'tif'); %save mask
    close all
    
end
