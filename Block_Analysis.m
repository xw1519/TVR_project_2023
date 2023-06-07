%% ANALYSIS BY EMG CHANNEL - FIXED TIME INTERVAL
nreps  = 5;                     % Number of repetitions of the trial
nChans = 14;                    % Number of channels EMG
Amp = [0, 85, 170, 255];        % Amplitudes considered
Optfreq = 120;                  % optimal frequency
tWindowInterval=[               % Window intervals to calculate delta EMG
   
4,5;
%{
    4,5;
    5,6;
    6,7;
    7,8;
%}
    6,8
    
%{
    10,11;
    11,12;
    12,13; 
    13,14;

    4,5;%stimulation starts
    13,14;
    14,15
%}
    ];

Path = '/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/Block_C_S3';

%initiate results
yData_High = [];
yData_Medium = [];
%yData_Low = [];
yData_Zero = [];

% Get a list of all files, including folders.
DirList  = dir(Path);
% Extract only the folders, not regular files.
DirList  = DirList([DirList.isdir]);  % Folders only
% Get rid of first two folders: dot and dot dot.
DirList = DirList(3:end);
% Warn user if there are no subfolders.
if isempty(DirList)
  message = sprintf('This folder does not contain any subfolders:\n%s', Path);
  uiwait(errordlg(message));
  return;
end
% Count the number of subfolders.
numberOfFolders = numel(DirList);
% Loop over all subfolders, processing each one.
for k = 1 : numberOfFolders
  thisDir = fullfile(Path, DirList(k).name);%new path
  fprintf('Processing folder %d of %d: %s\n', ...
    k, numberOfFolders, thisDir);
    
    fileList = dir(fullfile(thisDir, '/*_QUAT*.mat'));
    thisfile = fullfile(thisDir, fileList(1).name);
    load(thisfile);
    %Selected_channel = RecInfo.Calibration.ChannelsSelected(4); % char type
    %nchan = str2num(Selected_channel); %selected channel
    nchan = input('What are the manually identified optimal channel: '); %optimal channel

    T=[];
    % Calculate the gradient EMG for each channel and window
    for window = 2:length(tWindowInterval)
        fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
        T = [T; deltaEMG_3(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,thisDir)];
    end

    High_amp = T.Amplitude == 255;
    yData_High   = [yData_High; table2array(T(High_amp,2:end))];
    Medium_amp = T.Amplitude == 170;
    yData_Medium = [yData_Medium; table2array(T(Medium_amp,2:end))];
    %Low_amp = T.Amplitude == 85;
    %yData_Low    = [yData_Low; table2array(T(Low_amp,2:end))];
    %Zero_amp = T.Amplitude == 0;
    %yData_Zero   = [yData_Zero; table2array(T(Zero_amp,2:end))];
end
%disp(T1)
%% Gradient at the optimal channels across volunteers
figure
boxplot([yData_High,yData_Medium],{'High','Medium'})
ylabel('∆EMG'' [%]')
xlabel('Optimal Channel')
ylim([-4,10]);
grid
set(gca,'FontSize',24)
%% plotting - execute 1 line at one time 
GroupedData_0 = [yData_High yData_Medium yData_Low yData_Zero];
GroupedData_10 = [yData_High yData_Medium yData_Low yData_Zero];
GroupedData_20 = [yData_High yData_Medium yData_Low yData_Zero];
GroupedData_30 = [yData_High yData_Medium yData_Low yData_Zero];
%% boxplotgroup
x = {GroupedData_0, GroupedData_10, GroupedData_20, GroupedData_30};
figure('Position',[172 758 902 421])
boxplotGroup(x, 'Colors', lines(4))
title('A) BlockII Results - Amplitude / %MVC','FontName','FixedWidth')
set(gca,'FontSize',24, 'XTick', [])
colors = {[0 0 0] [0.00,0.34,0.74] [0.93,0.69,0.13] [0.49,0.18,0.56]};

grid on;
%% Gradient at the optimal channels across volunteers - color version

%window = '14<t<15';
%{
yData       = T.deltaEMG(strcmp(T.Period,window));
xData       = T.EMGChan(strcmp(T.Period,window));
ColorGroups = T.Amplitude(strcmp(T.Period,window));
%}
yData       = T.gradient;
%xData       = T.EMGChan;
ColorGroups = T.Amplitude;

F1 = figure(1);
boxchart(yData,'GroupByColor',ColorGroups);
ylabel('gradient ∆EMG'' [%]')
ylim([-1,4]);
grid
set(gca,'FontSize',24, 'XTick', [])
F1.Position=[172 758 902 421];
legend('zero','low','medium','high');
title('0%MVC Trial');
F1.Children(2).Children(4).BoxFaceColor=[0,0,0];
F1.Children(2).Children(3).BoxFaceColor=[0.00,0.34,0.74];
F1.Children(2).Children(2).BoxFaceColor=[0.93,0.69,0.13];
F1.Children(2).Children(1).BoxFaceColor=[0.49,0.18,0.56];