%% ANALYSIS BY EMG CHANNEL - FIXED TIME INTERVAL
nreps  = 5;                     % Number of repetitions of the trial
nChans = 14;                    % Number of channels EMG
Amp = [0, 85, 170, 255];        % Amplitudes considered
Optfreq = 120;                  % optimal frequency
tWindowInterval=[               % Window intervals to calculate delta EMG
    4,5;
    
    %4,5;
    %{
    5,6;
    6,7;
    %}
    6,8;
    %8,9;
    
    %{
    9,10;
    10,11;
    11,12;
    12,13; %before stimulation 
    
    
    13,14;
    14,15;
    
    
    15,16;
    16,17;
    17,18 %after stimulation
    %}
    ];

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Xinyao_force2N_0829_20']
};

T0_1=[];
% Calculate the delta EMG for each channel and window
for window = 2:length(tWindowInterval)
    %{
    for nchan = 1:nChans
        fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
        T0_1 = [T0_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
    end
    %}
    
    nchan = 4;
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_1 = [T0_1; deltaEMG_3(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
    
end

disp(T0_1)
%% Gradient over the channels

window = '14<t<15';
%{
yData       = T.deltaEMG(strcmp(T.Period,window));
xData       = T.EMGChan(strcmp(T.Period,window));g
ColorGroups = T.Amplitude(strcmp(T.Period,window));
%}
yData       = T0_1.deltaEMG;
xData       = T0_1.EMGChan;
ColorGroups = T0_1.Amplitude;

F1 = figure(1);
boxchart(xData, yData,'GroupByColor',ColorGroups);
ylabel('∆EMG'' [%]')
xlabel('Optimal Channel')
ylim([0,10]);
grid
set(gca,'FontSize',24)
F1.Position=[172 758 902 421];
F1.Children.XAxis.TickValues=[1:1:12];
legend('zero','low','medium','high');
F1.Children(2).Children(4).BoxFaceColor=[0,0,0];
F1.Children(2).Children(3).BoxFaceColor=[0.00,0.34,0.74];
F1.Children(2).Children(2).BoxFaceColor=[0.93,0.69,0.13];
F1.Children(2).Children(1).BoxFaceColor=[0.49,0.18,0.56];