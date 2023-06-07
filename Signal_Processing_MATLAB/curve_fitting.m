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

Chan = 11;                                       % define optimal channel

% AUX VARIABLES
cF = ones(1,length(Amps));                      % Counter per analysed trial
cF_2 = ones(1,length(Amps));

% PLOTTING VARIABLES
f = {};                                         % Variable to organise the plottings
% Color for plotting
color(1,:) = linspace(1,0,ntrials);
color(2,:) = zeros(1,ntrials);
color(3,:) = linspace(0,1,ntrials);

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
    fe = pF; % Assign peak frequency
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
    Fs_2 = find(Amps == VibAmp); 
    % Store the results on the defined position
    Amp(cF(Fs),Fs) = VibAmp;
    Amp(cF_2(Fs_2),Fs_2) = VibAmp;
    % Remove the time offset
    Time = RecInfo.Data.TimeStamp - ones(length(RecInfo.Data.TimeStamp),1) .* RecInfo.Data.TimeStamp(1);
    
    %explore individual trial variabilities
    MA_coef_num = 50; %smooth
    MA = ones(1,MA_coef_num)/MA_coef_num; %moving avaerage filter
    Data_filt = conv(Data, MA, 'same');
    
    M = movvar(Data, 5);%time window as 200ms - 500/20 No./Sec * 200/1000 = 5 
    M_enlarged = M * 100;
    f{Fs_2}=figure(Fs_2);
    minLength  = min (length(Time), length(Data));
    plot(Time(1:5:minLength),M(1:5:minLength)*100,'color',color(:,cF_2(Fs_2)));
    ylabel('Variance in a time window of 200ms');
    xlabel('Time [s]');
    %{
    % Consider data for time window 1s to 15s, calculate average
    Processed_data = zeros(1,150);
    count = 1;
    for timepoint = 1:0.1:15
        SamplesToAnalyse = find(Time > timepoint & Time < timepoint+0.1);
        Ave_data = Data(SamplesToAnalyse(1):SamplesToAnalyse(end));
        Processed_data (count) = mean(Ave_data);
        count = count + 1;
    end
    %}
    %{
    % plot original data
    f{Fs}=figure(Fs);
    plot(Time(1:minLength),Data(1:minLength).*100,'color', color(:,cF(Fs)))
    hold on
    plot(Time(1:minLength),Data_filt(1:minLength).*100,'Color', 'g','linewidth',2)
    hold on
    errorbar(Time(1:5:minLength),Data(1:5:minLength).*100,M_enlarged(1:5:minLength).*100,'LineStyle','none', 'color', color(:,cF(Fs)),'linewidth', 2);
    ylabel('Delta EMG Activity');
    xlabel('Time [s]');
    xlim([0,20]);
    ylim([0,60])
    set(gca,'fontsize',30);
    f{Fs}.Position=[175 630 1202 389];
    grid

    
    %}
    % fit linear for rate of activity increasing - not accurate 
    type lreval
    IncreaseToFit = find(Time > 4.5 & Time < 6);
    start_data = Data(IncreaseToFit(1));
    end_data = Data(IncreaseToFit(end));
    fit_time = Time(IncreaseToFit(1):IncreaseToFit(end));
    k = (end_data - start_data)/1.5;
    const = start_data - k*4.5;
    K_Set(cF(Fs),Fs) = k;

    % fit non-linear after activity rises
    type sseval
    SamplesToAnalyse = find(Time > 6 & Time < 9);
    analysed_data = Data(SamplesToAnalyse(1):SamplesToAnalyse(end));
    analysed_time = Time(SamplesToAnalyse(1):SamplesToAnalyse(end));
    ss_equation = @(x)sseval(x,analysed_time',analysed_data*100); 
    x0 = [1,1];
    bestx = fminsearch(ss_equation,x0);
    A_Set(cF(Fs),Fs) = bestx(1);
    B_Set(cF(Fs),Fs) = bestx(2);

    %check the fit quality
    A = bestx(1);
    b = bestx(2);
    f{Fs}=figure(Fs);
    hold on
    yfit = A*(1-exp(-analysed_time/b));
    increasefit = k*fit_time + const;
    minLength  = min (length(Time), length(Data));
    plot(Time(1:minLength),Data(1:minLength).*100,'color',color(:,cF(Fs)),'linewidth',2)
    hold on
    plot(analysed_time,yfit,'g');
    hold on
    plot(fit_time,increasefit.*100,'b');
    xlim([0,20]);
    ylim([0,30]);
    xlabel('Time [s]');
    ylabel('Normalized EMG siganl[au]')
    title('Data and Best Fitting Curve')
    legend('Data','Fitted Curve')
    hold off
    
    cF(Fs) = cF(Fs) + 1;

end

%Amplitude   = reshape(Amp,1,[])';
%A_value     = reshape(A_Set,1,[])';
%B_value     = reshape(B_Set,1,[])';

%T=table(Number,Amplitude,A_value,B_value);
%disp(T);