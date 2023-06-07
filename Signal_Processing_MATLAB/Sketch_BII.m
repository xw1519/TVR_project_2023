%Load mat dat a\
path=uigetdir();
ToScan = strcat(path,'/*_QUAT*.mat');%'in windows: \'
files = dir(ToScan);
SortedFiles = natsortfiles({files.name});    

% Initialize variables & arrays
Amps = [0, 85, 170, 255]; 
ntrials = 20;
tWindowInterval=[          % Window intervals to calculate delta Force
    4,5;% first evaluation
    %{
    4,5;
    5,6;
    6,7;
    7,8;
    8,9;
    9,10;
    10,11;
    11,12;
    12,13; %before stimulation 
    13,14;
    %} 
    14,15;% second evaluation
    %{
    15,16;
    16,17;
    17,18 %after stimulation
    %}
    18,19;% third evaluation
    ];

ControlConditions = length(Amps);               % Conditions to test
Force = zeros(ntrials,ControlConditions);       % Array to store the delta Force
WindowBlock = {};                               % Array to store the labels of the analysed window
log = zeros(ntrials,ControlConditions);         % Array to store all gradients
Amp = zeros(ntrials,ControlConditions);         % Array to store the labels per trial and amplitude

% AUX VARIABLES
cF = ones(1,length(Amps));                      % Counter per analysed trial
cF_2 = ones(1,20);    

% PLOTTING VARIABLES
f = {};                                         % Variable to organise the plottings

%Extract data from the files
for nfile=1:length(SortedFiles)
    % Load data from the *.mat file
    FileToLoad = fullfile(path,SortedFiles{nfile});
    load(FileToLoad);
    
    % Get trial stimulation amplitude, number, & force
    VibAmp = RecInfo.Experiment.Order(nfile);

    aux = RecInfo.Data.FilteredEMG'; 
    Data = aux(14,:);% chan14 EMG data is force
    % Remove the time offset
    Time = RecInfo.Data.TimeStamp - ones(length(RecInfo.Data.TimeStamp),1) .* RecInfo.Data.TimeStamp(1);

    % Calculate gradient for time window 5s to 15s
    SamplesToAnalyse = find(Time > 5 & Time < 15);
    y_Force = Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) ;
    X_time = Time(SamplesToAnalyse(1):SamplesToAnalyse(end));
    b = polyfit(X_time,y_Force, 1);
    b_gradient = b(1); % small gradient means stable

    % Get start part to analyse
    SamplesToAnalyse = find(Time > tWindowInterval(1,1) & Time < tWindowInterval(1,2));
    StartForce = rms( Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) );
        
    % Get end part to analyse
    SamplesToAnalyse = find(Time > tWindowInterval(2,1) & Time < tWindowInterval(2,2));
    EndForce = rms( Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) );
   
    % Calculate the delta EMG -  difference between beggining and end
    % of the trial
    difForce = EndForce - StartForce;

    % Extract Forces at three timepoints
    Force_1 = StartForce;
    Force_2 = EndForce;
    SamplesToAnalyse = find(Time > tWindowInterval(3,1) & Time < tWindowInterval(3,2));
    Force_3 = rms( Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) );

    % Find the index to store the results - amplitude 
    Fs = find(Amps == VibAmp); 
    % Store the results on the defined position
    Num(cF(Fs),Fs) = nfile;
    AForce(cF(Fs),Fs) = difForce;
    Amp(cF(Fs),Fs) = VibAmp;
    log(cF(Fs),Fs) = b_gradient;
    WindowBlock(cF(Fs),Fs) = {strcat(num2str(tWindowInterval(2,1)),'<t<',num2str(tWindowInterval(2,2)))};

    % Find the index to store the results - orders in trial
    Num_2(cF_2(nfile),nfile) = nfile;
    Amp_2(cF_2(nfile),nfile) = VibAmp;
    Force_start(cF_2(nfile),nfile) = Force_1;
    Force_middle(cF_2(nfile),nfile) = Force_2;
    Force_end(cF_2(nfile),nfile) = Force_3;

    %f{Fs}=figure(Fs);
    %hold on
    figure
    minLength  = min (length(Time), length(Data));
    plot(Time(1:minLength),Data(1:minLength))
    ylabel('Calibrated Stapping Force');
    xlabel('Time [s]');
    xlim([0,20]);
    ylim([0,3])
    set(gca,'fontsize',30);
    f{Fs}.Position=[175 630 1202 389];

end

%data arrangement - orders in trial version
Number      = reshape(Num_2,1,[])';
Amplitude   = reshape(Amp_2,1,[])';
Starting_Force    = reshape(Force_start,1,[])';
Middle_Force = reshape(Force_middle,1,[])';
Ending_Force = reshape(Force_end,1,[])';

T=table(Number,Amplitude,Starting_Force,Middle_Force,Ending_Force);
disp(T);

%%
% evaluate a single test
idx = T.Number < 2 & T.Number > 0;
T_single = T(idx,:);
Array = table2array(T_single);
A_adjusted = Array(1,3:5);
plot(A_adjusted);
ylim = ([0,3]);

%%
% evaluate in the whole block 
rows_medium = (T.Amplitude > 169 & T.Amplitude < 171);
T_MediumAmp = T.Middle_Force (rows_medium);
rows_high = (T.Amplitude > 171);
T_HighAmp = T.Middle_Force (rows_high);
figure
subplot(2,1,1);
plot(T_MediumAmp)
title('force changes with time - Medium Amplitude Trials')
xlabel('trial number')
ylabel('Calibrated Stapping Force')
subplot(2,1,2);
plot(T_HighAmp)
title('force changes with time - High Amplitude Trials')
xlabel('trial number')
ylabel('Calibrated Stapping Force')





