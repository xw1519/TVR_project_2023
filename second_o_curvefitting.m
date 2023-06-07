% Load mat dat a\
path=uigetdir();

% Get all data files
ToScan = strcat(path,'/*_QUAT*.mat');
files = dir(ToScan);
SortedFiles = natsortfiles({files.name});
% Load random file to detect the number of channels
load(fullfile(path,SortedFiles{1}))
nChannels = size(RecInfo.Data.rmsEMG,2);

% Initialize Variables
Amps = [0, 85, 170, 255]; 
ntrials = 20;
ControlConditions = length(Amps);
WindowBlock = {};                               % Array to store the labels of the analysed window
log = zeros(ntrials,ControlConditions);         % Array to store all gradients
Amp = zeros(ntrials,ControlConditions);         % Array to store the labels per trial and amplitude
AEMG = zeros(ntrials,ControlConditions);        % Array to store the delta EMG

Chan = 4;                                       % define optimal channel

% AUX VARIABLES
cF = ones(1,length(Amps));                      % Counter per analysed trial

% PLOTTING VARIABLES
f = {};                                         % Variable to organise the plottings
% Color for plotting
color(1,:) = linspace(1,0,ntrials);
color(2,:) = zeros(1,ntrials);
color(3,:) = linspace(0,1,ntrials);

% Initialize the cell arrays
data_1 = {};
data_2 = {};
data_3 = {};
data_4 = {};

%Extract data from the files
for nfile=1:length(SortedFiles)
    % Load the raw EMG from the BIN file
    BinName = strcat(SortedFiles{nfile}(1:end-3),'bin');               
    rawEMG = readBin_simple([nChannels,Inf],'int16',fullfile(path,BinName));
    % Load data from the *.mat file
    FileToLoad = fullfile(path,SortedFiles{nfile});
    load(FileToLoad);
    
    % Get trial stimulation amplitude
    VibAmp = RecInfo.Experiment.Order(nfile);

    % Raw data processed with notch filter
    % Detect noise peak frequency
    [pxx2,f2] = eH(rawEMG(Chan,:),1000,1,0);
    pM = max(abs(pxx2)); %magnitude
    pF = f2(abs(pxx2)==pM); %frequency

    rmsWindow = 500;
    BufSize = 40;
    windowCount = 1;
    rmsEMG = zeros(nChannels, round((length(rawEMG)/BufSize) - (rmsWindow/BufSize))) ;
    samplect = 1;

    % filter the noise at 120Hz
    Qfactor = 10;
    fe = 117; % Assign peak frequency
    nHarmonics = 2;
    fs = 1000; % system at 1k hz
    if fe>0
        frequencies = [1:nHarmonics].*fe;
        for freqToNotch = frequencies
            wo = freqToNotch/(fs/2);  
            bw = wo/Qfactor;
            [b,a] = iirnotch(wo,bw);
            rawEMG(Chan,:) = filtfilt(b,a,rawEMG(Chan,:));
        end
    end
    % filter the noise at 50Hz
    Qfactor = 10;
    fe = 49; % Assign peak frequency
    nHarmonics = 2;
    fs = 1000; % system at 1k hz
    if fe>0
        frequencies = [1:nHarmonics].*fe;
        for freqToNotch = frequencies
            wo = freqToNotch/(fs/2);  
            bw = wo/Qfactor;
            [b,a] = iirnotch(wo,bw);
            rawEMG(Chan,:) = filtfilt(b,a,rawEMG(Chan,:));
            %c=c+1;
        end
    end

    % Extract rms of the signal
    for n  = 1 : round((length(rawEMG)/BufSize) - (rmsWindow/BufSize))
        if n == 1 
            DataI = rawEMG(:, 1:rmsWindow);
        else
            DataI = rawEMG(:, (BufSize * (n-1) + 1):(BufSize * (n-1)) + rmsWindow);
        end
        % row, column
        rmsEMG(:,samplect) = rms(DataI');
        samplect = samplect + 1;
    end
    
    % Normalise data
    NrmsEMG = zeros(nChannels, round((length(rawEMG)/BufSize) - (rmsWindow/BufSize))) ;
    for nchan = 1 : nChannels
        NrmsEMG(nchan, :) = remap(rmsEMG(nchan, :),RecInfo.Calibration.EMVC(2,nchan),RecInfo.Calibration.EMVC(1,nchan),0,1);
    end
    
    % EMG channel to analyse
    aux = NrmsEMG; % offline analysis of the EMG
    %aux = RecInfo.Data.rmsNormEMG';

    Data = aux(Chan,:);% optimal chan only
    % Find the index to store the results - amplitude 
    Fs = find(Amps == VibAmp); 
    % Store the results on the defined position
    Amp(cF(Fs),Fs) = VibAmp;
    % Remove the time offset
    Time = RecInfo.Data.TimeStamp - ones(length(RecInfo.Data.TimeStamp),1) .* RecInfo.Data.TimeStamp(1);

    
    %This part is to plot variance quantification
    
    MA_coef_num = 50; 
    MA = ones(1,MA_coef_num)/MA_coef_num; %moving avaerage filter
    Data_filt = conv(Data, MA, 'same');
    %{
    M = movvar(Data, 5);%time window as 200ms - 500/20 No./Sec * 200/1000 = 5 
    M_enlarged = M * 100;
    f{Fs}=figure(Fs);
    minLength  = min (length(Time), length(Data));
    plot(Time(1:5:minLength),M(1:5:minLength)*100,'color',color(:,cF(Fs)));
    ylabel('Variance in a time window of 200ms');
    xlabel('Time [s]');
    %}

    % Consider data for time window 1s to 15s, calculate average
    %{
    Processed_data = zeros(1,150);
    count = 1;
    for timepoint = 0:0.1:15
        SamplesToAnalyse = find(Time > timepoint & Time < timepoint+0.1);
        Ave_data = Data(SamplesToAnalyse(1):SamplesToAnalyse(end));
        Processed_data (count) = mean(Ave_data);
        count = count + 1;
    end
    %}
   
    % store data
    switch Fs
        case 1
            data_1= [data_1; Data_filt];
        case 2
            data_2= [data_2; Data_filt];
        case 3
            data_3= [data_3; Data_filt];
        case 4
            data_4= [data_4; Data_filt];
    end

    % plot original data
    %{
    f{Fs}=figure(Fs);
    hold on
    minLength  = min (length(Time), length(Data_filt));
    plot(Time(1:minLength),Data_filt(1:minLength).*100,'color', color(:,cF(Fs)),'linewidth',2)
    ylabel('Norm. rms EMG [%]');
    xlabel('Time [s]');
    xlim([0,20]);
    ylim([0,30])
    set(gca,'fontsize',30);
    f{Fs}.Position=[175 630 1202 389];
    grid
    %}
    cF(Fs) = cF(Fs) + 1;
end

% equal length
for i = 1:length(data_1)
    data_1{i} = data_1{i}(1:450);
    data_4{i} = data_4{i}(1:450);
end

% Concatenate the truncated arrays along the third dimension
data_1_concat = cat(3, data_1{:});
data_4_concat = cat(3, data_4{:});

% mean curve and their difference
data_1_mean = mean(data_1_concat, 3);
data_4_mean = mean(data_4_concat, 3);
diff_data = zeros(size(data_1_mean));
for ii = 1:size(data_1_mean, 1)
    diff_data(ii,:) = data_4_mean(ii,:) - data_1_mean(ii,:);
end

% Second-Order Curve Fitting Plot
for ii = 1:size(data_1_mean, 1)
    figure(ii);
    plot(Time(1:450),data_1_mean(ii,:)*100, 'LineWidth', 2);
    xlim([3,10]);
    hold on;
    plot(Time(1:450),data_4_mean(ii,:)*100, 'LineWidth', 2);
    plot(Time(1:450),diff_data(ii,:)*100, 'LineWidth', 2);
    hold off;
    legend('Zero Amplitude', 'High Amplitude', 'Amplitude-Oriented Activity');
    ylabel('Norm. rms EMG [%]');
    xlabel('Time [s]');
    set(gca,'fontsize',16);
    f{Fs}.Position=[175 630 1202 389];
    grid
end