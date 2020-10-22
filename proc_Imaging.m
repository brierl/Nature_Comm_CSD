function [data_dot,gcamp6corr,WL,info,color_channels,gcamp6u]=proc_Imaging(rawFolder,filename,Date,gcamp_file,rs_gcamp6corr,rs_gcamp6u,rs_colorchannels)
% main processing script for mouse imaging data. Major steps include 1)
% subtract dark frame 2) downsample temporally 3) ratiometric correction if
% mice have gcamp 4) detrend, mean subtract, oximetry on hgb data 5) detrend, 
% mean subtract gcamp data if mice have gcamp 
% 6) get rid of NaN of Inf. If trial other than trial 1 is
% being processed, trial will have a CSD and greatly alter mean subtraction
% step. Therefore variables from baseline trial 1 run are loaded to perform
% mean subtraction.

% IN:
%   rawFolder: rawFolder/Date/ is directory containing raw data
%   filename: Assumes naming convention Date-Mouse-trials{i} for mouse 
%       imaging is used. i loops through run number. .tif file to process
%   Date: date mouse was imaged
%   gcamp_file: %flag for if mouse has gcamp, 1=yes, 0=no
%   rs_gcamp6corr: ratio corrected gcamp for hgb confound (Wright et al., 2017). 
%       if trial 1, empty variable, else, load gcamp6corr from trial 1
%   rs_gcamp6u: uncorrected gcamp for hgb confound.
%       if trial 1, empty variable, else, load gcamp6corr from trial 1
%   rs_colorchannels: individual LED channels for oximetry.
%       if trial 1, empty variable, else, load gcamp6corr from trial 1

% OUT:
%   data_dot=oxy(:,:,1,:) or deoxy(:,:,2,:) hgb. pixel x pixel x contrast x frame
%   gcamp6corr=ratio corrected gcamp for hgb confound (Wright et al., 2017). pixel x pixel x frame
%   WL=white light image. pixel x pixel
%   color_channels=individual LED channels for oximetry. pixel x pixel x LED# x frames
%   gcamp6u=uncorrected gcamp for hgb confound. pixel x pixel x frame
        
    disp('< < : ... L O A D I N G     Y O U R     D A T A ... : > >')
    
    %different settings if hgb imaging only vs hgb+gcamp
    if gcamp_file==1
        info.framerate=16.81;
    else
        info.framerate=29.41;
    end
   
    info.numled=4;
    info.freqout=1; %temporal resample will be to 1Hz
        
    %Get optical properties (Extinction/absorption coefficients, etc) and
    %optical instrument properties
    [op, E]=getop;
    
    %load raw data and resize so it can be later reshaped into 4 LED
    %channels
    [rawdata]=getoisdata([rawFolder Date '/' filename]);
    rawdata=double(rawdata);   
    [nVy, nVx, L]=size(rawdata);
    L2=L-rem(L,info.numled);
    rawdata=rawdata(:,:,1:L2);
    info.nVx=nVx; %number of pixels in x
    info.nVy=nVy; %number of pixels in y

    % 1) Subtract dark/background frame (ie no lights)
    dark=single(imread('Dark.tif')); 
    for z=1:L2
    	raw(:,:,z)=abs(rawdata(:,:,z)-dark);
    end
    clear rawdata
    
    %Reshape the data, adding a dimension for individiual LED channels
    rawdata=reshape(raw,info.nVy,info.nVx,info.numled,[]);
    clear raw
    %drop noisy first frame across all LED channels
    raw=rawdata(:,:,:,2:end); clear rawdata
     
    %%% Create white - light image:
    frameone=double(raw(:,:,:,2));
    WL(:,:,1)=frameone(:,:,2)/max(max(frameone(:,:,2)));
    WL(:,:,2)=frameone(:,:,3)/max(max(frameone(:,:,3)));
    WL(:,:,3)=frameone(:,:,4)/max(max(frameone(:,:,4)));

    % 2) resample temporally
    raw2=resampledata(raw,info.framerate,info.freqout,10^-5); 
   
    %Separate data into different channels. Rawdata is data from the green, yellow, and
    %red LEDs to be used to do oximetry. GCAMP6 is the raw fluorescence
    %emission channel. Green is the green channel on its own (used for the
    % 3) ratiometric fluroescence correction to remove hemodynamic confound)
    if gcamp_file==1 %if mice have gcamp
        rawdata=raw2(:,:,2:4,:);
        gcamp6=double(squeeze(raw2(:,:,1,:)));
        green=double(squeeze(raw2(:,:,2,:)));
        gcamp6c=gcamp6./green;
    else
        rawdata=raw2;
    end
    
    %prep to process pixels
    info.T1=size(rawdata,4); %number of frames
    info.numled=3; %number of LEDs for oximetry
    xsize=info.nVx; %number of pixels in x
    ysize=info.nVy; %number of pixels in y

    if size(rs_colorchannels,1)==0 %empty variable means this is resting state data, 1st trial of this mouse, can use same run for mean subtraction
        disp('< < : ... P R O C E S S I N G    H G B    R S ... : > >')
        
        for x=1:xsize %process one pixel at a time
            for y=1:ysize 
                % 4) "Process" the data--procPixel_rs used for detrending, mean subtraction, and oximetry   
                [data_dot(y,x,:,:), color_channels(y,x,:,:)]=procPixel_rs(squeeze(rawdata(y,x,:,:)), E);
            end
        end
        
        if gcamp_file==1 %if mice have gcamp
            disp('< < : ... P R O C E S S I N G    G C A M P   R S ... : > >')
            for x=1:xsize %process one pixel at a time
                for y=1:ysize        
                    
                    % 5) procPixel2_rs used for detrending and mean subtraction
                    gcamp6corr(y,x,:)=procPixel2_rs(squeeze(gcamp6c(y,x,:))'); %ratio corrected gcamp
                    gcamp6u(y,x,:)=procPixel2_rs(squeeze(gcamp6(y,x,:))'); %uncorrected gcamp
                
                end
            end
        end
        
    else %non-empty variable means this is a CSD run, must load 1st trial of this mouse to use for mean subtraction
        disp('< < : ... P R O C E S S I N G    H G B ... : > >')
        colors=rs_colorchannels; %baseline data per LED from 1st trial
        clear color_channels            
        
        for x=1:xsize %process one pixel at a time
            for y=1:ysize    
                % 4) "Process" the data--procPixel used for detrending, mean subtraction, and oximetry 
                [data_dot(y,x,:,:), color_channels(y,x,:,:)]=procPixel(squeeze(rawdata(y,x,:,:)), E, squeeze(colors(y,x,:,:)));
            end
        end
        
        if gcamp_file==1 %if mice have gcamp
            disp('< < : ... P R O C E S S I N G    G C A M P ... : > >')
            g=rs_gcamp6corr; %baseline ratio corrected gcamp data from 1st trial
            gu=rs_gcamp6u; %baseline uncorrected gcamp data from 1st trial
            for x=1:xsize %process one pixel at a time
                for y=1:ysize

                    % 5) procPixel2 used for detrending and mean subtraction
                    gcamp6corr(y,x,:)=procPixel2(squeeze(gcamp6c(y,x,:))', squeeze(g(y,x,:))); %ratio corrected gcamp
                    gcamp6u(y,x,:)=procPixel2(squeeze(gcamp6(y,x,:))', squeeze(gu(y,x,:))); %uncorrected gcamp
                    
                end
            end
        end
    end
    
    % 6) get rid of NaN or Inf in data
    data_dot(isnan(data_dot))=0;
    data_dot(isinf(data_dot))=0;
    
    if gcamp_file==1 %if mice have gcamp
    	gcamp6corr(isnan(gcamp6corr))=0;
        gcamp6corr(isinf(gcamp6corr))=0;    
        gcamp6u(isnan(gcamp6u))=0;
        gcamp6u(isinf(gcamp6u))=0;    
    else
        gcamp6corr=0;
        gcamp6u=0;
    end
    
end

%% getoisdata()
function [data]=getoisdata(fileall)
%read .tif stack (fileall) and assemble into matrix (data)

%get path, name, and type of file
[filepath,filename,filetype]=interppathstr(fileall);

if ~isempty(filetype) && ~(strcmp(filetype,'tif'))
    error('** script only supports the loading of .tif files **')
else
    data=readtiff([filepath,filename,'.tif']);
end

end

%% readtiff()
function [data]=readtiff(filename)
%read .tif stack (filename) and assemble into matrix (data)

info = imfinfo(filename);
numI = numel(info);
data=zeros(info(1).Width,info(1).Height,numI,'uint16');
fid=fopen(filename);

fseek(fid,info(1).Offset,'bof');
for k = 1:numI
    fseek(fid,[info(1,1).StripOffsets(1)-info(1).Offset],'cof');    
    tempdata=fread(fid,info(1).Width*info(1).Height,'uint16');
    data(:,:,k) = rot90((reshape(tempdata,info(1).Width,info(1).Height)),-1);
end

fclose(fid);

end

%% getop()
function [op, E, numled, led]=getop
% optical properties, spectroscopy for imaging system as described in
% Wright et al., 2017
% op: optical properties
% E: spectroscopy matrix
% numled: number of LEDs for oximetry
% led: spectra

[lambda1, Hb]=getHb; %get extinction coefficients (Hb) per wavelength (lambda1)
[led,lambda2]=getLED; %get LED spectra (led) per wavelength (lambda2)
   
op.HbT=76*10^-3; % uM concentration
op.sO2=0.71; % Oxygen saturation (%/100)
op.BV=0.1; % blood volume (%/100)

op.nin=1.4; % Internal Index of Refraction
op.nout=1; % External Index of Refraction
op.c=3e10/op.nin; % Speed of Light in the Medium
op.musp=10; % Reduced Scattering Coefficient

numled=size(led,2);

for n=1:numled                                                            
    
    % Interpolate from Spectrometer Wavelengths to Reference Wavelengths
    led{n}.ledpower=interp1(lambda2,led{n}.spectrum,lambda1,'pchip');
    
    % Normalize
    led{n}.ledpower=led{n}.ledpower/max(led{n}.ledpower);
    
    % Zero Out Noise
    led{n}.ledpower(led{n}.ledpower<0.01)=0;
    
    % Normalize
    led{n}.ledpower=led{n}.ledpower/sum(led{n}.ledpower);
    
    % Absorption Coeff.
    op.mua(n)=sum((Hb(:,1)*op.HbT*op.sO2+Hb(:,2)*op.HbT*(1-op.sO2)).*led{n}.ledpower);
    
    % Diffusion Coefficient
    op.gamma(n)=sqrt(op.c)/sqrt(3*(op.mua(n)+op.musp));
    op.dc(n)=1/(3*(op.mua(n)+op.musp));
    
    % Spectroscopy Matrix
    E(n,1)=sum(Hb(:,1).*led{n}.ledpower);
    E(n,2)=sum(Hb(:,2).*led{n}.ledpower);
    
    % Differential Pathlength Factors
    op.dpf(n)=(op.c/op.musp)*(1/(2*op.gamma(n)*sqrt(op.mua(n)*op.c)))*(1+(3/op.c)*op.mua(n)*op.gamma(n)^2);

end

end

function [lambda, Hb]=getHb
%get extinction coefficients (Hb) per wavelength (lambda)

data=dlmread('prahl_extinct_coef.txt'); %.txt file of extinction coefficients

lambda=data(:,1);
c=log(10)/10^3; %convert: (1) base-10 to base-e and (2) M^-1 to mM^-1
Hb=c*squeeze(data(:,2:3));

end

%% getLED()
function [led, lambda]=getLED
%get LED spectra (led) per wavelength (lambda)

led{1}.name='131029_Mightex_530nm_NoBPFilter'; %spectra files for LED used
led{2}.name='140801_ThorLabs_590nm_NoPol'; %spectra files for LED used
led{3}.name='140801_ThorLabs_625nm_NoPol'; %spectra files for LED used

numled=size(led,2);

%Read in LED spectra data from included text files
for n=1:numled
    
	fid=fopen([led{n}.name, '.txt']);
	temp=textscan(fid,'%f %f','headerlines',17);
	fclose(fid);
	lambda=temp{1};
	led{n}.spectrum=temp{2};
  
end

end

%% procPixel()
function [data_dot, data5]=procPixel(data,E,mean_pixels)
%if not first trial, CSDs might greatly effect mean subtraction, so pull
%from resting first trial

% IN: 
%   data: single pixel time trace
%   E: spectroscopy matrix
%   mean_pixels: single pixel time trace from trial 1

% OUT:
%   data_dot: hgb time trace for pixel (detrended and mean subtracted)
%   data5: detrended and mean subtracted time trace for pixel

for col=1:3 %loops through LEDs
    data4=detrend(data(col,:));
    data5(col,:)=data4-nanmean(mean_pixels(col,:)); 
    clear data4
end

%perform oximetry
data_dot=dotspect(data5,E(1:3,:));

end

function [data2]=procPixel2(data,mean_pixels)
%if not first trial, CSDs might greatly effect mean subtraction, so pull
%from resting first trial

% IN: 
%   data: single pixel time trace
%   mean_pixels: single pixel time trace from trial 1

% OUT:
%   data2: detrended and mean subtracted time trace for pixel

data4=detrend(data(:));
data2(:)=data4-nanmean(mean_pixels); 

end
%% procPixel_rs()
function [data_dot, data5]=procPixel_rs(data,E)
%if first trial

% IN: 
%   data: single pixel time trace
%   E: spectroscopy matrix

% OUT:
%   data_dot: hgb time trace for pixel (detrended and mean subtracted)
%   data5: detrended and mean subtracted time trace for pixel

for col=1:3 %loops through LEDs
    data4=detrend(data(col,:));
    data5(col,:)=data4-nanmean(data4); 
    clear data4
end

%perform oximetry
data_dot=dotspect(data5,E(1:3,:));

end

function [data2]=procPixel2_rs(data)
%if first trial

% IN: 
%   data: single pixel time trace

% OUT:
%   data2: detrended and mean subtracted time trace for pixel

data4=detrend(data(:));
data2(:)=data4-nanmean(data4); 

end
