% Load binary data
data=readBin_simple([16,Inf],'int16');

%% Load binary dat a\
path=uigetdir();
ToScan = strcat(path,'/*_QUAT*.bin');%'in windows: \'
files = dir(ToScan);
SortedFiles = natsortfiles({files.name});    

%%
clc
for nfile=1:length(SortedFiles)
    % Load the raw EMG from the BIN file
    BinName = strcat(SortedFiles{nfile}(1:end-3),'bin');               
    T_data{nfile} = readBin_simple([16,Inf],'int16',fullfile(path,BinName));
    disp(BinName)
end
%% Plot all channels
trial = 3;
data = T_data{trial};
ChansToPlot = [4];
Fs = 1000;
PlotChs(data, ChansToPlot ,Fs)

%% Filtering
Fs = 1000; %sampling f
Fc = 20;%cutoff f
chan = 4;
order = 9; %steep, great for cutoff

% Low pass filter of the signal
[b,a] = butter(order,Fc/(Fs/2),'low');%normalized cutoff f
dataFilt = filtfilt(b,a,data(chan,:));

% Plotting
figure(1)
plot(data(chan,:),'color',[1,0,0],'linewidth',1)
hold on
plot(dataFilt,'color',[0,0,0],'linewidth',1)
plot(data(14,:) + 200,'color',[0,0,0],'linewidth',1) 

% Check in power spectrum if signal was filtered
figure(2)
eH(dataFilt,1000,1)

%% Spectrum analysis
timeToAnalyse = [10,16] .*1000;
Channels = 4;
Fs = 1000;
Smoothing = 0.5;

%data = data_LS;

% Plot the spectrum
figure(2)
eH(data(Channels,timeToAnalyse(1):timeToAnalyse(2)),Fs,Smoothing);

labels = {''};
for i = 1:length(Channels)
    labels{i} = strcat('Chan',num2str(Channels(i)));
end
legend(labels)
set(gca,'Fontsize',30)

