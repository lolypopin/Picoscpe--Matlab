%%%% NONEED TO CHANGE
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
autoTriggerMicroSeconds = 1000000; %us

[status.setSimpleTrigger] = invoke(ps6000aDeviceObj, 'ps6000aSetSimpleTrigger', enable, source, threshold, direction,...
    delay, autoTriggerMicroSeconds);

disp('Simple Trigger set')


%% Get Fastest Timebase 4

enabledChannelFlags= ps6000aEnumInfo.enPicoChannelFlags.PICO_CHANNEL_A_FLAGS + ps6000aEnumInfo.enPicoChannelFlags.PICO_CHANNEL_B_FLAGS;
pTimebase = libpointer('uint32Ptr',0);
pTimeInterval = libpointer('doublePtr',0);

[status.getMinimumTimebaseStateless] = invoke(ps6000aDeviceObj, 'ps6000aGetMinimumTimebaseStateless', enabledChannelFlags,...
    pTimebase, pTimeInterval, resolution);

timebase = pTimebase.Value;
timeInterval = pTimeInterval.Value;


    
    
    

%datasetnumber=1000000;  %%%DATASET NUMBER 1,000,000
e=zeros(100,1);
t2=zeros(100,1);
ecount=1;
disp("ecount");
datasetnumber = linspace(2500000,250000000, 100);
for loop=1:100
    disp("loop");
    disp(loop);
    nSegments = 10;
    nMaxSamples = datasetnumber(loop); %1,000,000
    pnMaxSamples = libpointer('uint64Ptr', nMaxSamples);
    [status.memorySegments] = invoke(ps6000aDeviceObj, 'ps6000aMemorySegments', nSegments, pnMaxSamples);

    %% Set number of samples to be collected
    
    numPreTriggerSamples = 0;
    numPostTriggerSamples = nMaxSamples/10; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    totalSamples = numPreTriggerSamples + numPostTriggerSamples;
    
    %% Set number of captures
    
    [status.setNoOfCaptures] = invoke(ps6000aDeviceObj, 'ps6000aSetNoOfCaptures', nSegments);
    
    bufferAMax = zeros(totalSamples, 1, 'int16');
    bufferBMax = zeros(totalSamples, 1, 'int16');
    
    for i=(1:10)
        pBufferAMax(i) =libpointer('int16Ptr', bufferAMax);
        pBufferBMax(i) =libpointer('int16Ptr', bufferBMax);
    end
    
    dataType = ps6000aEnumInfo.enPicoDataType.PICO_INT16_T;
    downSampleRatioMode = ps6000aEnumInfo.enPicoRatioMode.PICO_RATIO_MODE_AVERAGE;
    actionA = bitor(ps6000aEnumInfo.enPicoAction.PICO_CLEAR_ALL, ps6000aEnumInfo.enPicoAction.PICO_ADD);
    actionB = ps6000aEnumInfo.enPicoAction.PICO_ADD;
    
    %% Create Buffers  8

    [status.setBufferA.zero] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(1), ...
        totalSamples, dataType, 0, downSampleRatioMode, actionA);
    [status.setBufferB.zero] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(1), ...
        totalSamples, dataType, 0, downSampleRatioMode, actionB);
    
    [status.setBufferA.one] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(2), ...
        totalSamples, dataType, 1, downSampleRatioMode, actionB);
    [status.setBufferB.one] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(2), ...
        totalSamples, dataType, 1, downSampleRatioMode, actionB);
    
    [status.setBufferA.two] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(3), ...
        totalSamples, dataType, 2, downSampleRatioMode, actionB);
    [status.setBufferB.two] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(3), ...
        totalSamples, dataType, 2, downSampleRatioMode, actionB);
    
    [status.setBufferA.three] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(4), ...
        totalSamples, dataType, 3, downSampleRatioMode, actionB);
    [status.setBufferB.three] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(4), ...
        totalSamples, dataType, 3, downSampleRatioMode, actionB);
    
    [status.setBufferA.four] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(5), ...
        totalSamples, dataType, 4, downSampleRatioMode, actionB);
    [status.setBufferB.four] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(5), ...
        totalSamples, dataType, 4, downSampleRatioMode, actionB);
    
    [status.setBufferA.five] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(6), ...
        totalSamples, dataType, 5, downSampleRatioMode, actionB);
    [status.setBufferB.five] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(6), ...
        totalSamples, dataType, 5, downSampleRatioMode, actionB);
    
    [status.setBufferA.six] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(7), ...
        totalSamples, dataType, 6, downSampleRatioMode, actionB);
    [status.setBufferB.six] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(7), ...
        totalSamples, dataType, 6, downSampleRatioMode, actionB);
    
    [status.setBufferA.seven] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(8), ...
        totalSamples, dataType, 7, downSampleRatioMode, actionB);
    [status.setBufferB.seven] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(8), ...
        totalSamples, dataType, 7, downSampleRatioMode, actionB);
    
    [status.setBufferA.eight] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(9), ...
        totalSamples, dataType, 8, downSampleRatioMode, actionB);
    [status.setBufferB.eight] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(9), ...
        totalSamples, dataType, 8, downSampleRatioMode, actionB);
    
    [status.setBufferA.nine] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelA, pBufferAMax(10), ...
        totalSamples, dataType, 9, downSampleRatioMode, actionB);
    [status.setBufferB.nine] = invoke(ps6000aDeviceObj, 'ps6000aSetDataBuffer', channelB, pBufferBMax(10), ...
        totalSamples, dataType, 9, downSampleRatioMode, actionB);
    
    %% Run Block Capture
    
    pTimeIndisposedMs = libpointer('doublePtr',0);
    segmentIndex = 0;
    
    disp('Collection starting...')
    
    [status.runBlock] = invoke(ps6000aDeviceObj, 'ps6000aRunBlock', numPreTriggerSamples, numPostTriggerSamples,...
        timebase, pTimeIndisposedMs, segmentIndex); 
    
    pReady = libpointer('int16Ptr',0);
    
    while pReady.Value == 0
        [status.IsReady] = invoke(ps6000aDeviceObj,'ps6000aIsReady',pReady);
    end
    
    disp('Collection finished')
    
    %% Retrieve Data  9

    startIndex = 0;
    pSamplesCollected = libpointer('uint64Ptr',totalSamples);
    downSampleRatio = 1;
    segmentIndex = 0;
    pOverflow = libpointer('int16Ptr',zeros(10,1));
    fromSegmentIndex = 0;
    toSegmentIndex = 9;
    
    [status.getValuesBulk] = invoke(ps6000aDeviceObj, 'ps6000aGetValuesBulk', startIndex,...
        pSamplesCollected, fromSegmentIndex, toSegmentIndex, downSampleRatio, downSampleRatioMode, pOverflow);
    
    samplesCollected = pSamplesCollected.Value;
    
    disp('Data Retrieved')
    
    %% Convert Data from ADC counts to mV
    
    BufferAMax={};
    BufferBMax={};
    for i=(1:10)
        BufferAMax{i} = pBufferAMax(i).Value;
        BufferBMax{i} = pBufferBMax(i).Value;
    end
    
    pMinValue = libpointer('int16Ptr',0);
    pMaxValue = libpointer('int16Ptr',0);
    
    [status.getAdcLimits] = invoke(ps6000aDeviceObj, 'ps6000aGetAdcLimits', resolution, pMinValue, pMaxValue);
    
    minValue = pMinValue.Value;
    maxValue = pMaxValue.Value;
    
    voltageRange = 5000; %mV
    
    bufferAMax={};
    bufferBMax={};
    
    for i=(1:10)
        bufferAMax{i} = adc2mv(BufferAMax{i},voltageRange,double(maxValue));
        bufferBMax{i} = adc2mv(BufferBMax{i},voltageRange,double(maxValue));
    end
    
    disp('Data converted to mV')
    
    BufferAMax={};
    BufferBMax={};
    for i=(1:10)
        BufferAMax{i} = pBufferAMax(i).Value;
        BufferBMax{i} = pBufferBMax(i).Value;
    end
    
    pMinValue = libpointer('int16Ptr',0);
    pMaxValue = libpointer('int16Ptr',0);
    
    [status.getAdcLimits] = invoke(ps6000aDeviceObj, 'ps6000aGetAdcLimits', resolution, pMinValue, pMaxValue);
    
    minValue = pMinValue.Value;
    maxValue = pMaxValue.Value;
    
    voltageRange = 5000; %mV
    
    bufferAMax={};
    bufferBMax={};
    
    for i=(1:10)
        bufferAMax{i} = adc2mv(BufferAMax{i},voltageRange,double(maxValue));
        bufferBMax{i} = adc2mv(BufferBMax{i},voltageRange,double(maxValue));
    end
    
    disp('Data converted to mV')
    
    %% Plot Collected Data
    
    maxTime = (double(samplesCollected) * timeInterval);
    time = linspace(0,maxTime,samplesCollected);
    
    time3= linspace(0,maxTime*4,samplesCollected*4);
    totalBuffers = [];
    n_cal=4;
    for i = 1:n_cal
        totalBuffers = vertcat(totalBuffers,bufferBMax{1,i});
    end
    %% Calculate the time of the gap
    
    % Prerequisite of variables
    %
    %-bufferB is the chunks of data collected.

    %Set the number of C (sample_mean) as 1/10 of a chunk
    sample_mean=samplesCollected/10.0;
    disp(sample_mean);
    
    % INPUT Configuration
    Vpp=1.00; % Enter Vpp of input signal
    Frequency=50;% Enter the frequency of input signal 
    
    % create an array for the phase shifted data
    phaseshifted_bufferBMax=bufferBMax;
    phaseshifted_bufferBMax_mean={};
    shift_count=0;
    shifted=0;

    % Phase shift No.2 for 2 chunks
    for j=1:2
        for i=1:sample_mean
            % if reset occured in C; Phase shift
            if bufferBMax{1,j}(1)*bufferBMax{1,j}(i) <0 && bufferBMax{1,j}(1)>0.4*Vpp
                phaseshifted_bufferBMax{1,j}(i)=bufferBMax{1,j}(i)+Vpp*1000;
            end
            phaseshifted_bufferBMax_mean{1,j}(i,1)=phaseshifted_bufferBMax{1,j}(i);
        end
    end

    disp( phaseshifted_bufferBMax{1,1}(1) );
    disp(length(phaseshifted_bufferBMax_mean{1,1}));

    % Create phase shifted data for plotting
    shiftedBuffers = [];
    shiftedBuffers1 = [];
    for i = 1:2
        shiftedBuffers = vertcat(shiftedBuffers, phaseshifted_bufferBMax{1,i});
        % The length of C
        shiftedBuffers1 = vertcat(shiftedBuffers1, phaseshifted_bufferBMax_mean{1, i});
    end

    disp(length(shiftedBuffers));
    % Create a time vector for 2 Cs
    tt1=2 * (double(sample_mean) * timeInterval);
    time2= linspace(0,tt1,sample_mean*2);
    
    %plot C1 and C2 
    figure(4)
    hold on
    plot(time2, shiftedBuffers1);
    ylabel('Voltage (mV)');
    xlabel('Time (s)');
    hold off
    
    mean=zeros(2, 1);
    % Calculate the mean of C1 and C2
    for j=1:2
        for i=1:sample_mean
            mean(j)=double(mean(j)+phaseshifted_bufferBMax_mean{1, j}(i));
        end
    end
    mean(1)=mean(1)/double(sample_mean);
    mean(2)=mean(2)/double(sample_mean);
    
    % when mean(1)>mean2 -> phase shift
    if mean(1)>mean(2)
        disp("1")
        disp(mean(2))
        mean(2)=mean(2)+Vpp*1000;
        disp(mean(2))
    end
    
    % Calculate D
    t2_1=mean(2)-mean(1);
    disp("t2_1")
    disp(t2_1)

    %Calculate time gap t2 using the Step 2
    t2_1=(1/double(Frequency))*(t2_1/(Vpp*1000));
    disp("t2_1- calculated")
    disp(t2_1)
    disp(timeInterval*double(samplesCollected))

    % Calculate the t1 using Step1
    t1=timeInterval*double(samplesCollected);

    % Calculate the efficiency using Step3
    efficiency=double(t1)/double(t2_1);
    disp(efficiency)
    e(ecount)=efficiency;
    t2(ecount)=t2_1;
    ecount=ecount+1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Disconnect scope

disconnect(ps6000aDeviceObj);

%%
delete(ps6000aDeviceObj);