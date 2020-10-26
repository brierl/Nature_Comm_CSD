function fcOIS_GCaMP()

% run fcOIS imaging system for ST seconds

ST=300; %seconds of recording per run
AO={'ao0','ao5','ao2','ao3','ao4','ao7','ao1'}; %NI daq board analog outputs CCD, LED1, LED2, LED3, LED4, EEG, Stimbox
Sr=10000; %sampling rate, default 10000
total_framerate=68; %framerate per LED x num of LEDs, default 68
Stim='Off'; %'On' or 'Off'
BlockLength=60; %length (s) of one stim block if Stim='On' or resting state block if Stim='Off'
BL1=5; %length (s) of baseline recording pre stim if Stim='On', doesn't matter if Stim='Off'
Dur=10; %duration of stim (s) if Stim='On', doesn't matter if Stim='Off'

%% AO control
outV=5;     % output trigger height (volts)

s = daq.createSession('ni');
s.Rate = Sr;
Ar = s.Rate;

addAnalogOutputChannel(s, 'Dev1', AO{1}, 'Voltage');%CCD
addAnalogOutputChannel(s, 'Dev1', AO{2}, 'Voltage');%High powered LED 470nm
addAnalogOutputChannel(s, 'Dev1', AO{3}, 'Voltage');%LED2 530nm
addAnalogOutputChannel(s, 'Dev1', AO{4}, 'Voltage');%LED3 590nm
addAnalogOutputChannel(s, 'Dev1', AO{5}, 'Voltage');%LED4 625nm
addAnalogOutputChannel(s, 'Dev1', AO{6}, 'Voltage');%EEG
addAnalogOutputChannel(s, 'Dev1', AO{7}, 'Voltage');%Stim

%% LED Settings
TExp=1/total_framerate; %exposure time
FR=1/TExp; %framerate
SpF=floor(Ar/FR); %sample rate/framerate

LED_SpF=[SpF SpF SpF SpF]; %Camera exposure length
LED_time=[146 20 13 3]; %LED duty cycle

LED_num=size(LED_time,2); %num of LEDs
scalfac=sum(LED_SpF)/(LED_num*SpF);

clear DataVec2

LEDind=1;
LEDTrigs=zeros(SpF, LED_num);

for n=1:LED_num %loop through LEDs, assign triggers for different LED exposure lengths for one camera exposure for each LED
    
    if LED_time(n)>=LED_SpF(n)
        error(['LED_', num2str(n), ' duration greater than camera exposure time.']);
    end
    
    startInd=round((SpF-round(LED_time(n)))./2); %start LED trigger
    endInd=startInd+round(LED_time(n)); %end LED trigger
    
    LEDTrigs(startInd:endInd, LEDind)=outV; %assign voltage during trigger
    LEDind=LEDind+1;
    clear DataVec_LED
end

DataVec2=zeros(SpF*LED_num, LED_num+1); %camera exposure time*num LEDs for one BGYR frame. Add another indice in the second dim for camera trigger

for n=1:LED_num
       
    startInd=SpF*(n-1)+1; %assign a camera trigger for each LED trigger
    DataVec2(startInd, 1)=outV; %fill in camera trigger for DataVec2
    endInd=startInd+SpF-1; %find indice for each end of camera exposure
    DataVec2(startInd:endInd, n+1)=LEDTrigs(:,n); %fill DataVec2 with LED triggers
    
end
   
%% Stimulus Settings
switch Stim
        
    case 'On'           % Use for stimulating in a block design
        BL2=BlockLength-(BL1+Dur); % Time between stim off and start of next block (s)
 
    case 'Off'              % No Stimuli
        BL1=BlockLength;
        Dur=0;
        BL2=0;
 
end

DataVec2(:,end+1)=outV; %add trigger for starting EEG recording
DataVecStim=DataVec2;
DataVecStim(:,end+1)=outV; %add trigger if Stim='On'
DataVecNoStim=DataVec2;
DataVecNoStim(:,end+1)=0; %no trigger if Stim='Off'

%% Concatenate Everything Together
DataVecOff1=repmat(DataVecNoStim,round(BL1*Sr/(LED_num*SpF*scalfac)),1); %repmat without stim trigger if Stim='Off' or Stim='On' for length Bl1*Sr/(LED_num*SpF*scalfac
DataVecOn=repmat(DataVecStim,round(Dur*Sr/(LED_num*SpF*scalfac)),1); %repmat with stim trigger if Stim='On' for length Dur*Sr/(LED_num*SpF*scalfac
DataVecOff2=repmat(DataVecNoStim,round(BL2*Sr/(LED_num*SpF*scalfac)),1); %repmat without stim trigger if Stim='On' for length Bl2*Sr/(LED_num*SpF*scalfac

DataVecTotal=cat(1, DataVecOff1, DataVecOn, DataVecOff2); %concatanate stim blocks together, or resting state blocks together
Repeats=single(round(ST/(BL1+Dur+BL2))); %find how many stim/resting state blocks fit in ST
DataVecTotal=repmat(DataVecTotal,Repeats,1); %repmat for Repeats

%% Begin Scanning

disp(['Set camera to record ',num2str(ST*FR),' frames.'])
disp('Press ENTER when ready to start flashing sequence.')
pause
disp('Playing encoding sequence.')
queueOutputData(s, DataVecTotal);
startForeground(s);

end
