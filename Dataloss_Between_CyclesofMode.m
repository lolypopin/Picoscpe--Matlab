%% Clear Command Window and Close Figures

clc;
close all;

%% Load Configuration Information

[ps6000aStructs, ps6000aEnumInfo]=PS6000aSetConfig();

%% Device Connection

% Check if an Instrument session using the device object 'ps6000DeviceObj'
% is still open, and if so, disconnect if the User chooses 'Yes' when prompted.
if (exist('ps6000aDeviceObj', 'var') && ps6000aDeviceObj.isvalid && strcmp(ps6000aDeviceObj.status, 'open'))
    
    openDevice = questionDialog(['Device object ps6000aDeviceObj has an open connection. ' ...
        'Do you wish to close the connection and continue?'], ...
        'Device Object Connection Open');
    
    if (openDevice == PicoConstants.TRUE)
        
        % Close connection to device
        disconnect(ps6000aDeviceObj);
        delete(ps6000aDeviceObj);
        
    else
    
        % Exit script if User selects 'No'
        return;
        
    end
    
end

%% Create a device object. 
% The serial number can be specified as a second input parameter.

ps6000aDeviceObj = icdevice('picotech_ps6000a_generic.mdd','');

%% Connect scope

connect(ps6000aDeviceObj)

%% Set Device Resolution

resolution = ps6000aEnumInfo.enPicoDeviceResolution.PICO_DR_10BIT;

[status.setResolution] = invoke(ps6000aDeviceObj, 'ps6000aSetDeviceResolution', resolution);

disp('Device Resolution set to 10 bits')

%% Enable Channel A + B
% Disable  other channels
for i = (0:7)
    try
        [status.setChannelOff] = invoke(ps6000aDeviceObj, 'ps6000aSetChannelOff', i);
    catch
        
    end 
end

for j = (128:1:131)
    try
        [status.turnDigitalPortOff] = invoke(ps6000aDeviceObj, 'ps6000aDigitalPortOff',j);
    catch
    end
end

% Enable channels A + B with +-5 V range with DC coupling and full bandwidth

channelA = ps6000aEnumInfo.enPicoChannel.PICO_CHANNEL_A;
channelB = ps6000aEnumInfo.enPicoChannel.PICO_CHANNEL_B;
couplingDC = ps6000aEnumInfo.enPicoCoupling.PICO_DC_50OHM;
range = ps6000aEnumInfo.enPicoConnectProbeRange.PICO_X1_PROBE_5V;
bandwidth = ps6000aEnumInfo.enPicoBandwidthLimiter.PICO_BW_FULL;


[status.setChannelOn.A] = invoke(ps6000aDeviceObj, 'ps6000aSetChannelOn', channelA, couplingDC, range, 0, bandwidth);
[status.setChannelOn.B] = invoke(ps6000aDeviceObj, 'ps6000aSetChannelOn', channelB, couplingDC, range, 0, bandwidth);

disp('Channels A and B set')

%% Set Simple Trigger

enable = 1;
source = channelA;
threshold = 1000; %mV
direction = ps6000aEnumInfo.enPicoThresholdDirection.PICO_RISING;
delay = 0;
autoTriggerMicroSeconds = 1; %us

[status.setSimpleTrigger] = invoke(ps6000aDeviceObj, 'ps6000aSetSimpleTrigger', enable, source, threshold, direction,...
    delay, autoTriggerMicroSeconds);

disp('Simple Trigger set')


%% Get Fastest Timebase

enabledChannelFlags= ps6000aEnumInfo.enPicoChannelFlags.PICO_CHANNEL_A_FLAGS + ps6000aEnumInfo.enPicoChannelFlags.PICO_CHANNEL_B_FLAGS;
pTimebase = libpointer('uint32Ptr',0);
pTimeInterval = libpointer('doublePtr',0);

[status.getMinimumTimebaseStateless] = invoke(ps6000aDeviceObj, 'ps6000aGetMinimumTimebaseStateless', enabledChannelFlags,...
    pTimebase, pTimeInterval, resolution);

timebase = pTimebase.Value;
timeInterval = pTimeInterval.Value;

e_array=zeros(100,1);%efficiency array
t2_array=zeros(100, 1);
datasetnumber=zeros(100,1);
datasetnumber(1)=1000000;
for i =2:100
    datasetnumber(i)=datasetnumber(i-1)+50000;
end
for loop=1:100
    %% Set number of samples to be collected
    
    numPreTriggerSamples = 0;
    numPostTriggerSamples = datasetnumber(loop); % 100,000
    totalSamples = numPreTriggerSamples + numPostTriggerSamples;
    
    %% Create Buffers
        bufferAMax = zeros(totalSamples, 1, 'int16');
        bufferBMax = zeros(totalSamples, 1, 'int16');
        bufferAMin = zeros(totalSamples, 1, 'int16');
        bufferBMin = zeros(totalSamples, 1, 'int16');
        pBufferAMax =libpointer('int16Ptr', bufferAMax);
        pBufferBMax =libpointer('int16Ptr', bufferBMax);
        pBufferAMin =libpointer('int16Ptr', bufferAMin);
        pBufferBMin =libpointer('int16Ptr', bufferBMin);
    
        dataType = ps6000aEnumInfo.enPicoDataType.PICO_INT16_T;
        waveform = 0;
        downSampleRatioMode = ps6000aEnumInfo.enPicoRatioMode.PICO_RATIO_MODE_AVERAGE;
        actionA = bitor(ps6000aEnumInfo.enPicoAction.PICO_CLEAR_ALL, ps6000aEnumInfo.enPicoAction.PICO_ADD);
        actionB = ps6000aEnumInfo.enPicoAction.PICO_ADD;
    
        [status.setBufferA] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffers', channelA, pBufferAMax, ...
            pBufferAMin, totalSamples, dataType, waveform, downSampleRatioMode, actionA);
        [status.setBufferB] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffers', channelB, pBufferBMax, ...
            pBufferBMin, totalSamples, dataType, waveform, downSampleRatioMode, actionB);
    tic
    
    for p=1:2
        %% Run Block Capture (5)
        pTimeIndisposedMs = libpointer('doublePtr',0);
        segmentIndex = 0;
    
        disp('Collection starting...')
        toc
        [status.runBlock] = invoke(ps6000aDeviceObj, 'ps6000aRunBlock', numPreTriggerSamples, numPostTriggerSamples,...
            timebase, pTimeIndisposedMs, segmentIndex); 
    
        pReady = libpointer('int16Ptr',0);
        % BlockIsReady (6)
        while pReady.Value == 0
            [status.IsReady] = invoke(ps6000aDeviceObj,'ps6000aIsReady',pReady);
        end
        disp('Collection finished')
    
        %% Retrieve Data
        startIndex = 0;
        pSamplesCollected = libpointer('uint64Ptr',totalSamples);
        downSampleRatio = 1;
        segmentIndex = 0;
        pOverflow = libpointer('int16Ptr',0);
        % getValues (8) --send to usb
        [status.getValues] = invoke(ps6000aDeviceObj,'ps6000aGetValues', startIndex,...
            pSamplesCollected, downSampleRatio, downSampleRatioMode, segmentIndex, pOverflow);
    
        samplesCollected = pSamplesCollected.Value;
        tic
        disp('Data Retrieved')
        
    
        %% Convert Data from ADC counts to mV
        if p==1
            bufferAMax1 = pBufferAMax.Value;
            bufferAMin1 = pBufferAMin.Value;
            bufferBMax1 = pBufferBMax.Value;
            bufferBMin1 = pBufferBMin.Value;
        
            pMinValue1 = libpointer('int16Ptr',0);
            pMaxValue1 = libpointer('int16Ptr',0);
        
            [status.getAdcLimits] = invoke(ps6000aDeviceObj, 'ps6000aGetAdcLimits', resolution, pMinValue1, pMaxValue1);
        
            minValue1 = pMinValue1.Value;
            maxValue1 = pMaxValue1.Value;
        end
        if p==2
            bufferAMax2 = pBufferAMax.Value;
            bufferAMin2 = pBufferAMin.Value;
            bufferBMax2 = pBufferBMax.Value;
            bufferBMin2 = pBufferBMin.Value;
        
            pMinValue2 = libpointer('int16Ptr',0);
            pMaxValue2 = libpointer('int16Ptr',0);
        
            [status.getAdcLimits] = invoke(ps6000aDeviceObj, 'ps6000aGetAdcLimits', resolution, pMinValue2, pMaxValue2);
        
            minValue2 = pMinValue2.Value;
            maxValue2 = pMaxValue2.Value;
        end
    end
    
        voltageRange = 5000; %mV
        
        bufferAMax1 = adc2mv(bufferAMax1,voltageRange,double(maxValue1));
        bufferBMax1 = adc2mv(bufferBMax1,voltageRange,double(maxValue1));
    
        bufferAMax2 = adc2mv(bufferAMax2,voltageRange,double(maxValue2));
        bufferBMax2 = adc2mv(bufferBMax2,voltageRange,double(maxValue2));
        disp('Data converted to mV')
        
    
    %% Plot Collected Data
    listItemObjOutA = [bufferAMax1; bufferAMax2];
    listItemObjOutB = [bufferBMax1; bufferBMax2];
    
    maxTime = (double(samplesCollected) * timeInterval);
    time = linspace(0,2*maxTime,2*samplesCollected);
    
    
    figure(1)
    plot(time,listItemObjOutA);
    hold on
    plot(time,listItemObjOutB);
    ylabel('Voltage (mV)');
    xlabel('Time (s)');
    hold off
    
    sample_mean=samplesCollected/10.0;
    disp(sample_mean);
    mean=zeros(3, 1);
    disp(mean(1));
    Vpp=1; % Enter Vpp of Input
    F=2; %Enter Frequency of Input
    
    phaseshifted_bufferBMax1=bufferBMax1;
    phaseshifted_bufferBMax2=bufferBMax2;
    disp("#########")
    disp(length(bufferBMax1))
    disp(samplesCollected)
    %phase shift
    for i=0:samplesCollected-1
        if bufferBMax1(1)*bufferBMax1(i+1) <0 && bufferBMax1(1)>0.4*Vpp
            phaseshifted_bufferBMax1(i)=phaseshifted_bufferBMax1(i)+Vpp*1000;
        end
        if bufferBMax2(1)*bufferBMax2(i+1) <0 && bufferBMax2(1)>0.4*Vpp
            phaseshifted_bufferBMax2(i)=phaseshifted_bufferBMax2(i)+Vpp*1000;
        end
    end
    
    % calculate mean
    for i=1:sample_mean
        mean(1)=mean(1)+phaseshifted_bufferBMax1(i);
        mean(2)=mean(2)+phaseshifted_bufferBMax2(i);
    end
    mean(1)=mean(1)/double(sample_mean);
    mean(2)=mean(2)/double(sample_mean);
    %disp mean
    disp(mean(1));
    disp(mean(2));
    
    if mean(1) > mean(2)
        mean(2)= mean(2)+Vpp*1000;
        disp("!")%
    end
    
    %calc D
    D=mean(2)-mean(1);
    t2=(1/double(F))*(D/(Vpp*1000)); % T/num of sample
    t2_array(loop)=t2;
    t1=maxTime;
    efficiency=t1/t2;
    e_array(loop)=efficiency;
    disp("Efficieny")%
    disp(efficiency)
    disp(loop)%

end
%% Disconnect scope

disconnect(ps6000aDeviceObj);

%%

delete(ps6000aDeviceObj);