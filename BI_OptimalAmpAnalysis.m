%% ANALYSIS BY EMG CHANNEL - FIXED TIME INTERVAL
nreps  = 5;                     % Number of repetitions of the trial
nChans = 14;                    % Number of channels EMG
Amp = [0, 85, 170, 255];        % Amplitudes considered
Optfreq = 120;                  % optimal frequency
tWindowInterval=[               % Window intervals to calculate delta EMG
    %4,5;
    
    4,5;
    %{
    5,6;
    6,7;
    7,8;
    %}
    8,9;

    9,10;
    %{
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
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/0323_binghuan_1+2']
};

T0_1=[];
% Calculate the delta EMG for each channel and window
for window = 2:length(tWindowInterval)
    for nchan = 1:nChans
        fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
        T0_1 = [T0_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
    end
end

disp(T0_1)

%% Extract data only for optimal channels - FIXED TIME INTERVAL
nreps  = 5;                     % Number of repetitions of the trial
nchan = 14;                      % Optimal channel is 8
Amp = [0, 85, 170, 255];        % Amplitudes considered
tWindowInterval=[               % Window intervals to calculate delta EMG
    %{
    4,5;
    4,5;
    5,6;
    %}
    6,7;
    7,8;
    8,9;
    9,10;
    10,11;
    11,12;
    12,13; %before stimulation 
    13,14;
    %{
    14,15;
    15,16;
    16,17;
    17,18 %after stimulation
    %}
    ];

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pat0808']
};

T0_1=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_1 = [T0_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

disp(T0_1)

nchan = 9;

T0_1_monitorised=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_1_monitorised = [T0_1_monitorised; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol0810_00']
};

nchan = 8; 

T0_2=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_2 = [T0_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Zimu0812_00']
};

nchan = 8; 

T0_3=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_3 = [T0_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 9;

T0_3_monitorised=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_3_monitorised = [T0_3_monitorised; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pat0808_10']
};

nchan = 8; 

T1_1=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_1 = [T1_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 9;

T1_1_monitorised=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_1_monitorised = [T1_1_monitorised; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol0810_10']
};

nchan = 8; 

T1_2=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_2 = [T1_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Zimu0812_10']
};

nchan = 8; 

T1_3=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_3 = [T1_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 9;

T1_3_monitorised=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_3_monitorised = [T1_3_monitorised; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,9,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol0810_20']
};

nchan = 8; 

T2_2=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_2 = [T2_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Zimu0812_20']
};

nchan = 8; 

T2_3=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_3 = [T2_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 9;

T2_3_monitorised=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_3_monitorised = [T2_3_monitorised; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol0810_30']
};

nchan = 8; 

T3_2=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T3_2 = [T3_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Zimu0812_30']
};

nchan = 8; 

T3_3=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T3_3 = [T3_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 9; 

T3_3_monitorised=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T3_3_monitorised = [T3_3_monitorised; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

%plot 'difference between deltaMVC' figure
record_00 = [T0_1.(6)(6:20); T0_1.(6)(26:40); T0_2.(6)(6:20); T0_2.(6)(26:40); T0_3.(6)(6:20); T0_3.(6)(26:40)];
record_10 = [T1_1.(6)(6:20); T1_1.(6)(26:40); T1_2.(6)(6:20); T1_2.(6)(26:40); T1_3.(6)(6:20); T1_3.(6)(26:40)];
record_20 = [T2_2.(6)(6:20); T2_2.(6)(26:40); T2_3.(6)(6:20); T2_3.(6)(26:40)];
record_30 = [T3_2.(6)(6:20); T3_2.(6)(26:40); T3_3.(6)(6:20); T3_3.(6)(26:40)];

yData = [record_00;record_10;record_20;record_30];
group_yData = [ ones(size(record_00));
            2* ones(size(record_10));
            3* ones(size(record_20));
            4* ones(size(record_30))];


record_00_monitorised = [T0_1_monitorised.(6)(6:20); T0_1_monitorised.(6)(26:40); T0_2.(6)(6:20); T0_2.(6)(26:40); T0_3_monitorised.(6)(6:20); T0_3_monitorised.(6)(26:40)];
record_10_monitorised = [T1_1_monitorised.(6)(6:20); T1_1_monitorised.(6)(26:40); T1_2.(6)(6:20); T1_2.(6)(26:40); T1_3_monitorised.(6)(6:20); T1_3_monitorised.(6)(26:40)];
record_20_monitorised = [T2_2.(6)(6:20); T2_2.(6)(26:40); T2_3_monitorised.(6)(6:20); T2_3_monitorised.(6)(26:40)];
record_30_monitorised = [T3_2.(6)(6:20); T3_2.(6)(26:40); T3_3_monitorised.(6)(6:20); T3_3_monitorised.(6)(26:40)];

yData_2 = [record_00_monitorised; record_10_monitorised; record_20_monitorised; record_30_monitorised];
group_yData_2 = [ ones(size(record_00_monitorised));
            2* ones(size(record_10_monitorised));
            3* ones(size(record_20_monitorised));
            4* ones(size(record_30_monitorised))];

%% %MVC Comparison

F1 = figure(1);%optimal
boxplot(yData, group_yData);
ylabel('gradient ∆EMG [%]')
grid
set(gca,'FontSize',24, 'XTickLabel', {'0% MVC','10% MVC','20% MVC','30% MVC'})
F1.Position=[172 758 902 421];
title('%MVC at Optimal Channel (Zero Amplitude Trials Removed)','FontSize',24, 'FontWeight', 'normal', 'FontName', 'Calibri')
%legend('optimal','monitorised');

F2 = figure(2);%monitorised
boxplot(yData_2, group_yData_2);
ylabel('gradient ∆EMG [%]')
grid
set(gca,'FontSize',24, 'XTickLabel', {'0% MVC','10% MVC','20% MVC','30% MVC'})
F2.Position=[172 758 902 421];
%legend('optimal','monitorised');

%% High VS Median VS Zero Amplitude Comparison 
%case (2+2+3+3)

%high = '255';
record_00_high = [T0_1.(6)(16:20); T0_1.(6)(36:40); T0_2.(6)(16:20); T0_2.(6)(36:40); T0_3.(6)(16:20); T0_3.(6)(36:40)];
record_10_high = [T1_1.(6)(16:20); T1_1.(6)(36:40); T1_2.(6)(16:20); T1_2.(6)(36:40); T1_3.(6)(16:20); T1_3.(6)(36:40)];
record_20_high = [T2_2.(6)(16:20); T2_2.(6)(36:40); T2_3.(6)(16:20); T2_3.(6)(36:40)];
record_30_high = [T3_2.(6)(16:20); T3_2.(6)(36:40); T3_3.(6)(16:20); T3_3.(6)(36:40)];

yData_high = [record_00_high; record_10_high; record_20_high; record_30_high];
group_yData_high = [ ones(size(record_00_high));
            2* ones(size(record_10_high));
            3* ones(size(record_20_high));
            4* ones(size(record_30_high))];

%median = '170';
record_00_median = [T0_1.(6)(11:15); T0_1.(6)(31:35); T0_2.(6)(11:15); T0_2.(6)(31:35); T0_3.(6)(11:15); T0_3.(6)(31:35)];
record_10_median = [T1_1.(6)(11:15); T1_1.(6)(31:35); T1_2.(6)(11:15); T1_2.(6)(31:35); T1_3.(6)(11:15); T1_3.(6)(31:35)];
record_20_median = [T2_2.(6)(11:15); T2_2.(6)(31:35); T2_3.(6)(11:15); T2_3.(6)(31:35)]; 
record_30_median = [T3_2.(6)(11:15); T3_2.(6)(31:35); T3_3.(6)(11:15); T3_3.(6)(31:35)];

yData_median = [record_00_median;record_10_median;record_20_median;record_30_median];
group_yData_median = [ ones(size(record_00_median));
            2* ones(size(record_10_median));
            3* ones(size(record_20_median));
            4* ones(size(record_30_median))];

%zero = '0';
record_00_zero = [T0_1.(6)(1:5); T0_1.(6)(1:5); T0_2.(6)(1:5); T0_2.(6)(1:5); T0_3.(6)(1:5); T0_3.(6)(1:5)];
record_10_zero = [T1_1.(6)(1:5); T1_1.(6)(1:5); T1_2.(6)(1:5); T1_2.(6)(1:5); T1_3.(6)(1:5); T1_3.(6)(1:5)];
record_20_zero = [T2_2.(6)(1:5); T2_2.(6)(1:5); T2_3.(6)(1:5); T2_3.(6)(1:5)];
record_30_zero = [T3_2.(6)(1:5); T3_2.(6)(1:5); T3_3.(6)(1:5); T3_3.(6)(1:5)];

yData_zero = [record_00_zero; record_10_zero; record_20_zero; record_30_zero];
group_yData_zero = [ ones(size(record_00_zero));
            2* ones(size(record_10_zero));
            3* ones(size(record_20_zero));
            4* ones(size(record_30_zero))];

figure
subplot(2,3,1)
boxplot(yData_high, group_yData_high);
ylabel('∆EMG gradient [%]')
ylim([-2,4])
title('High Amplitude at Optimal Channel','FontSize',12, 'FontWeight', 'normal', 'FontName', 'Calibri')
grid
set(gca,'FontSize',10, 'XTickLabel', {'0% MVC','10% MVC','20% MVC','30% MVC'})
%F1.Position=[172 758 902 421];
%legend('optimal','monitorised');

subplot(2,3,2)
boxplot(yData_median, group_yData_median);
ylabel('∆EMG gradient [%]')
ylim([-2,4])
title('Medium Amplitude at Optimal Channel','FontSize',12, 'FontWeight', 'normal', 'FontName', 'Calibri')
grid
set(gca,'FontSize',10, 'XTickLabel', {'0% MVC','10% MVC','20% MVC','30% MVC'})
%F2.Position=[172 758 902 421];

subplot(2,3,3)
boxplot(yData_zero, group_yData_zero);
ylabel('∆EMG gradient [%]')
title('Zero Amplitude at Optimal Channel','FontSize',12, 'FontWeight', 'normal', 'FontName', 'Calibri')
ylim([-2,4])
grid
set(gca,'FontSize',10, 'XTickLabel', {'0% MVC','10% MVC','20% MVC','30% MVC'})
%F2.Position=[172 758 902 421];

%For moneterised channels
%high = '255';
record_00_high_monitorised = [T0_1_monitorised.(6)(16:20); T0_1_monitorised.(6)(36:40); T0_2.(6)(16:20); T0_2.(6)(36:40); T0_3_monitorised.(6)(16:20); T0_3_monitorised.(6)(36:40)];
record_10_high_monitorised = [T1_1_monitorised.(6)(16:20); T1_1_monitorised.(6)(36:40); T1_2.(6)(16:20); T1_2.(6)(36:40); T1_3_monitorised.(6)(16:20); T1_3_monitorised.(6)(36:40)];
record_20_high_monitorised = [T2_2.(6)(16:20); T2_2.(6)(36:40); T2_3_monitorised.(6)(16:20); T2_3_monitorised.(6)(36:40)];
record_30_high_monitorised = [T3_2.(6)(16:20); T3_2.(6)(36:40); T3_3_monitorised.(6)(16:20); T3_3_monitorised.(6)(36:40)];

yData_high_monitorised = [record_00_high_monitorised; record_10_high_monitorised; record_20_high_monitorised; record_30_high_monitorised];
group_yData_high_monitorised = [ ones(size(record_00_high_monitorised));
            2* ones(size(record_10_high_monitorised));
            3* ones(size(record_20_high_monitorised));
            4* ones(size(record_30_high_monitorised))];

%median = '170';
record_00_median_monitorised = [T0_1_monitorised.(6)(11:15); T0_1_monitorised.(6)(31:35); T0_2.(6)(11:15); T0_2.(6)(31:35); T0_3_monitorised.(6)(11:15); T0_3_monitorised.(6)(31:35)];
record_10_median_monitorised = [T1_1_monitorised.(6)(11:15); T1_1_monitorised.(6)(31:35); T1_2.(6)(11:15); T1_2.(6)(31:35); T1_3_monitorised.(6)(11:15); T1_3_monitorised.(6)(31:35)];
record_20_median_monitorised = [T2_2.(6)(11:15); T2_2.(6)(31:35); T2_3.(6)(11:15); T2_3_monitorised.(6)(31:35)]; 
record_30_median_monitorised = [T3_2.(6)(11:15); T3_2.(6)(31:35); T3_3.(6)(11:15); T3_3_monitorised.(6)(31:35)];

yData_median_monitorised = [record_00_median_monitorised;record_10_median_monitorised;record_20_median_monitorised;record_30_median_monitorised];
group_yData_median_monitorised = [ ones(size(record_00_median_monitorised));
            2* ones(size(record_10_median_monitorised));
            3* ones(size(record_20_median_monitorised));
            4* ones(size(record_30_median_monitorised))];

%high = '0';
record_00_zero_monitorised = [T0_1_monitorised.(6)(1:5); T0_1_monitorised.(6)(1:5); T0_2.(6)(1:5); T0_2.(6)(1:5); T0_3_monitorised.(6)(1:5); T0_3_monitorised.(6)(1:5)];
record_10_zero_monitorised = [T1_1_monitorised.(6)(1:5); T1_1_monitorised.(6)(1:5); T1_2.(6)(1:5); T1_2.(6)(1:5); T1_3_monitorised.(6)(1:5); T1_3_monitorised.(6)(1:5)];
record_20_zero_monitorised = [T2_2.(6)(1:5); T2_2.(6)(1:5); T2_3_monitorised.(6)(1:5); T2_3_monitorised.(6)(1:5)];
record_30_zero_monitorised = [T3_2.(6)(1:5); T3_2.(6)(1:5); T3_3_monitorised.(6)(1:5); T3_3_monitorised.(6)(1:5)];

yData_zero_monitorised = [record_00_zero_monitorised; record_10_zero_monitorised; record_20_zero_monitorised; record_30_zero_monitorised];
group_yData_zero_monitorised = [ ones(size(record_00_zero_monitorised));
            2* ones(size(record_10_zero_monitorised));
            3* ones(size(record_20_zero_monitorised));
            4* ones(size(record_30_zero_monitorised))];

subplot(2,3,4)
boxplot(yData_high_monitorised, group_yData_high_monitorised);
ylabel('∆EMG gradient [%]')
grid
set(gca,'FontSize',10, 'XTickLabel', {'0% MVC','10% MVC','20% MVC','30% MVC'})
%F1.Position=[172 758 902 421];
%legend('optimal','monitorised');

subplot(2,3,5)
boxplot(yData_median_monitorised, group_yData_median_monitorised);
ylabel('∆EMG gradient [%]')
grid
set(gca,'FontSize',10, 'XTickLabel', {'0% MVC','10% MVC','20% MVC','30% MVC'})
%F2.Position=[172 758 902 421];

subplot(2,3,6)
boxplot(yData_zero_monitorised, group_yData_zero_monitorised);
ylabel('∆EMG gradient [%]')
grid
set(gca,'FontSize',10, 'XTickLabel', {'0% MVC','10% MVC','20% MVC','30% MVC'})
%F2.Position=[172 758 902 421];
%legend('high','median');

%{
%case only considering two volunteers
record_00_high = [T0_2.gradient(16:20,6); T0_2.gradient(36:40,6); T0_3.gradient(16:20,6); T0_3.gradient(36:40,6)];
record_10_high = [T1_2.gradient(16:20,6); T1_2.gradient(36:40,6); T1_3.gradient(16:20,6); T1_3.gradient(36:40,6)];
record_20_high = [T2_2.gradient(16:20,6); T2_2.gradient(36:40,6); T2_3.gradient(16:20,6); T2_3.gradient(36:40,6)];
record_30_high = [T3_2.gradient(16:20,6); T3_2.gradient(36:40,6); T3_3.gradient(16:20,6); T3_3.gradient(36:40,6)];
yData_high = [record_00_high; record_10_high; record_20_high; record_30_high];

record_00_median = [T0_1(11:15,6); T0_1(31:35,6); T0_2(11:15,6); T0_2(31:35,6); T0_3(11:15,6); T0_3(31:35,6)];
record_10_median = [T1_1(11:15,6); T1_1(31:35,6); T1_2(11:15,6); T1_2(31:35,6); T1_3(11:15,6); T1_3(31:35,6)];
record_20_median = [T2_2(11:15,6); T2_2(31:35,6); T2_3(11:15,6); T2_3(31:35,6)]; 
record_30_median = [T3_2(11:15,6); T3_2(31:35,6); T3_3(11:15,6); T3_3(31:35,6)];
yData_median = [record_00_median;record_10_median;record_20_median;record_30_median];

F1 = figure(1);
boxplot(yData_high);
ylabel('∆EMG gradient [%]')
grid
set(gca,'FontSize',24, 'XTickLabel', {'0% MVC','10% MVC','20% MVC','30% MVC'})
F1.Position=[172 758 902 421];
%legend('optimal','monitorised');
%}

%% Extract data only for optimal channels - FIXED TIME INTERVAL - Contact Area 4 subjects
nreps  = 5;                     % Number of repetitions of the trial
nchan = 9;                      % Optimal channel is 9
Amp = [0, 85, 170, 255];        % Amplitudes considered
tWindowInterval=[               % Window intervals to calculate delta EMG
    4,5;
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
    %}
    13,14;
    14,15;
    %{
    15,16;
    16,17;
    17,18 %after stimulation
    %}
    ];

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pilot_15cm_0831_00']
};

T0_m_1=[];% 0%MVC, 1st vol, medium 
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_m_1 = [T0_m_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

disp(T0_m_1)

nchan = 9;

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pilot_15cm_0831_10']
};
T1_m_1=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_m_1 = [T1_m_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 10; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pilot_15cm_0831_20']
};

T2_m_1=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_m_1 = [T2_m_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 8; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Xinyao_0829_small_00']
};

T0_s_1=[]; % 0%MVC, 1st vol, small 
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_s_1 = [T0_s_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Xinyao_0829_small_10']
};

T1_s_1=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_s_1 = [T1_s_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Xinyao_0829_small_20']
};

T2_s_1=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_s_1 = [T2_s_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 10;

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Fei_medium_0830_00']
};

T0_m_2=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_m_2 = [T0_m_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 9; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Fei_medium_0830_10']
};

T1_m_2=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_m_2 = [T1_m_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 11; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Fei_medium_0830_20']
};

T2_m_2=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_m_2 = [T2_m_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 7;

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Fei_small_0830_00']
};

T0_s_2=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_s_2 = [T0_s_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,9,Path{1})];
end

nchan = 10;

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Fei_small_0830_10']
};

T1_s_2=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_s_2 = [T1_s_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Fei_small_0830_20']
};

T2_s_2=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_s_2 = [T2_s_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 8;

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pat_0905_00_2']
};

T0_m_3=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_m_3 = [T0_m_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 11;

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pat_0905_10']
};

T1_m_3=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_m_3 = [T1_m_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 10;

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pat_0905_20']
};

T2_m_3=[];
% Calculate the delta EMG for the optimal channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_m_3 = [T2_m_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 7; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pat_0905_00_9cm']
};

T0_s_3=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_s_3 = [T0_s_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 10; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pat_0905_10_9cm']
};

T1_s_3=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_s_3 = [T1_s_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 7; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Pat_0905_20_9cm']
};

T2_s_3=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_s_3 = [T2_s_3; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 11; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol_0905_00_15']
};

T0_m_4=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_m_4 = [T0_m_4; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 8; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol_0905_10_15']
};

T1_m_4=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_m_4 = [T1_m_4; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol_0905_20_15']
};

T2_m_4=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_m_4 = [T2_m_4; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 5; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol_0905_00_9']
};

T0_s_4=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_s_4 = [T0_s_4; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 8; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol_0905_10_9']
};

T1_s_4=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T1_s_4 = [T1_s_4; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 6; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Vol_0905_20_9']
};

T2_s_4=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T2_s_4 = [T2_s_4; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

%extract data from medium contact area - 15cm
record_00_medium_high = [T0_m_1.(6)(16:20); T0_m_1.(6)(36:40); T0_m_2.(6)(16:20); T0_m_2.(6)(36:40);T0_m_3.(6)(16:20); T0_m_3.(6)(36:40);T0_m_4.(6)(16:20); T0_m_4.(6)(36:40)];
record_00_medium_medium = [T0_m_1.(6)(11:15); T0_m_1.(6)(31:35); T0_m_2.(6)(11:15); T0_m_2.(6)(31:35);T0_m_3.(6)(11:15); T0_m_3.(6)(31:35);T0_m_4.(6)(11:15); T0_m_4.(6)(31:35)];
record_00_medium_zero = [T0_m_1.(6)(1:5); T0_m_1.(6)(21:25); T0_m_2.(6)(1:5); T0_m_2.(6)(21:25);T0_m_3.(6)(1:5); T0_m_3.(6)(21:25);T0_m_4.(6)(1:5); T0_m_4.(6)(21:25)];

yData_00_medium = [record_00_medium_high;record_00_medium_medium;record_00_medium_zero];

record_10_medium_high = [T1_m_1.(6)(16:20); T1_m_1.(6)(36:40); T1_m_2.(6)(16:20); T1_m_2.(6)(36:40);T1_m_3.(6)(16:20); T1_m_3.(6)(36:40);T1_m_4.(6)(16:20); T1_m_4.(6)(36:40)];
record_10_medium_medium = [T1_m_1.(6)(11:15); T1_m_1.(6)(31:35); T1_m_2.(6)(11:15); T1_m_2.(6)(31:35);T1_m_3.(6)(11:15); T1_m_3.(6)(31:35);T1_m_4.(6)(11:15); T1_m_4.(6)(31:35)];
record_10_medium_zero = [T1_m_1.(6)(1:5); T1_m_1.(6)(21:25); T1_m_2.(6)(1:5); T1_m_2.(6)(21:25);T1_m_3.(6)(1:5); T1_m_3.(6)(21:25);T1_m_4.(6)(1:5); T1_m_4.(6)(21:25)];

yData_10_medium = [record_10_medium_high;record_10_medium_medium;record_10_medium_zero];

record_20_medium_high = [T2_m_1.(6)(16:20); T2_m_1.(6)(36:40); T2_m_2.(6)(16:20); T2_m_2.(6)(36:40);T2_m_3.(6)(16:20); T2_m_3.(6)(36:40);T2_m_4.(6)(16:20); T2_m_4.(6)(36:40)];
record_20_medium_medium = [T2_m_1.(6)(11:15); T2_m_1.(6)(31:35); T2_m_2.(6)(11:15); T2_m_2.(6)(31:35);T2_m_3.(6)(11:15); T2_m_3.(6)(31:35);T2_m_4.(6)(11:15); T2_m_4.(6)(31:35)];
record_20_medium_zero = [T2_m_1.(6)(1:5); T2_m_1.(6)(21:25); T2_m_2.(6)(1:5); T2_m_2.(6)(21:25);T2_m_3.(6)(1:5); T2_m_3.(6)(21:25);T2_m_4.(6)(1:5); T2_m_4.(6)(21:25)];

yData_20_medium = [record_20_medium_high;record_20_medium_medium;record_20_medium_zero];

%extract data from small contact area - 9cm
record_00_small_high = [T0_s_1.(6)(16:20); T0_s_1.(6)(36:40); T0_s_2.(6)(16:20); T0_s_2.(6)(36:40);T0_s_3.(6)(16:20); T0_s_3.(6)(36:40);T0_s_4.(6)(16:20); T0_s_4.(6)(36:40)];
record_00_small_medium = [T0_s_1.(6)(11:15); T0_s_1.(6)(31:35); T0_s_2.(6)(11:15); T0_s_2.(6)(31:35);T0_s_3.(6)(11:15); T0_s_3.(6)(31:35);T0_s_4.(6)(11:15); T0_s_4.(6)(31:35)];
record_00_small_zero = [T0_s_1.(6)(1:5); T0_s_1.(6)(21:25); T0_s_2.(6)(1:5); T0_s_2.(6)(21:25);T0_s_3.(6)(1:5); T0_s_3.(6)(21:25);T0_s_4.(6)(1:5); T0_s_4.(6)(21:25)];

yData_00_small = [record_00_small_high;record_00_small_medium;record_00_small_zero];

record_10_small_high = [T1_s_1.(6)(16:20); T1_s_1.(6)(36:40); T1_s_2.(6)(16:20); T1_s_2.(6)(36:40);T1_s_3.(6)(16:20); T1_s_3.(6)(36:40);T1_s_4.(6)(16:20); T1_s_4.(6)(36:40)];
record_10_small_medium = [T1_s_1.(6)(11:15); T1_s_1.(6)(31:35); T1_s_2.(6)(11:15); T1_s_2.(6)(31:35);T1_s_3.(6)(11:15); T1_s_3.(6)(31:35);T1_s_4.(6)(11:15); T1_s_4.(6)(31:35)];
record_10_small_zero = [T1_s_1.(6)(1:5); T1_s_1.(6)(21:25); T1_s_2.(6)(1:5); T1_s_2.(6)(21:25);T1_s_3.(6)(1:5); T1_s_3.(6)(21:25);T1_s_4.(6)(1:5); T1_s_4.(6)(21:25)];

yData_10_small = [record_10_small_high;record_10_small_medium;record_10_small_zero];

record_20_small_high = [T2_s_1.(6)(16:20); T2_s_1.(6)(36:40); T2_s_2.(6)(16:20); T2_s_2.(6)(36:40);T2_s_3.(6)(16:20); T2_s_3.(6)(36:40);T2_s_4.(6)(16:20); T2_s_4.(6)(36:40)];
record_20_small_medium = [T2_s_1.(6)(11:15); T2_s_1.(6)(31:35); T2_s_2.(6)(11:15); T2_s_2.(6)(31:35);T2_s_3.(6)(11:15); T2_s_3.(6)(31:35);T2_s_4.(6)(11:15); T2_s_4.(6)(31:35)];
record_20_small_zero = [T2_s_1.(6)(1:5); T2_s_1.(6)(21:25); T2_s_2.(6)(1:5); T2_s_2.(6)(21:25);T2_s_3.(6)(1:5); T2_s_3.(6)(21:25);T2_s_4.(6)(1:5); T2_s_4.(6)(21:25)];

yData_20_small = [record_20_small_high;record_20_small_medium;record_20_small_zero];

%% Block II Contact Area

% a big table T
T_medium_00 = [T0_m_1;T0_m_2;T0_m_3;T0_m_4];
T_medium_10 = [T1_m_1;T1_m_2;T1_m_3;T1_m_4];
T_medium_20 = [T2_m_1;T2_m_2;T2_m_3;T2_m_4];
T_small_00  = [T0_s_1;T0_s_2;T0_s_3;T0_s_4];
T_small_10  = [T1_s_1;T1_s_2;T1_s_3;T1_s_4];
T_small_20  = [T2_s_1;T2_s_2;T2_s_3;T2_s_4];

% plot medium contact area - 15cm
grid on
subplot(1,2,1);
boxchart(T_medium_10.gradient,'GroupByColor',T_medium_10.Amplitude);
ylabel('Gradient ∆EMG'' [%]', 'FontSize',16)
xlabel('15cm Diameter','FontSize',16)
ylim([-6,10]);
grid on
subplot(1,2,2);
boxchart(T_small_10.gradient,'GroupByColor',T_small_10.Amplitude);
ylabel('Gradient ∆EMG'' [%]', 'FontSize',16)
xlabel('9cm Diameter','FontSize',16)
ylim([-6,10]);
grid on
legend('zero','low','medium','high');
%{

subplot(2,3,3);
boxchart(T_medium_20.gradient,'GroupByColor',T_medium_20.Amplitude);
ylabel('Gradient ∆EMG'' [%]', 'FontSize',16)
xlabel('20% MVC','FontSize',16)
ylim([-1,2.5]);
grid on
legend('zero','low','medium','high');

% plot small contact area - 9cm
grid on
subplot(2,3,4);
boxchart(T_small_00.gradient,'GroupByColor',T_small_00.Amplitude);
ylabel('Gradient ∆EMG'' [%]', 'FontSize',16)
xlabel('0% MVC','FontSize',16)
ylim([-1,2.5]);
grid on
subplot(2,3,5);
boxchart(T_small_10.gradient,'GroupByColor',T_small_10.Amplitude);
ylabel('Gradient ∆EMG'' [%]', 'FontSize',16)
xlabel('10% MVC','FontSize',16)
ylim([-1,2.5]);
grid on
subplot(2,3,6);
boxchart(T_small_20.gradient,'GroupByColor',T_small_20.Amplitude);
ylabel('Gradient ∆EMG'' [%]', 'FontSize',16)
xlabel('20% MVC','FontSize',16)
ylim([-1,2.5]);
grid on
legend('zero','low','medium','high');
%}

%{
F1.Position=[172 758 902 421];
F1.Children.XAxis.TickValues=[1:1:12];
legend('zero','low','medium','high');
F1.Children(2).Children(4).BoxFaceColor=[0,0,0];
F1.Children(2).Children(3).BoxFaceColor=[0.00,0.34,0.74];
F1.Children(2).Children(2).BoxFaceColor=[0.93,0.69,0.13];
F1.Children(2).Children(1).BoxFaceColor=[0.49,0.18,0.56];
%}
%title('Medium Contact Area - 15cm Diameter Round Shape','FontSize',24, 'FontWeight', 'normal', 'FontName', 'Calibri')
%set(gca,'FontSize',24)

%%
%strapping force verification trials
nchan = 10; 

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Xinyao_force2N_0829_20']
};

T_2N=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T_2N = [T_2N; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end

nchan = 9;
Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/Xinyao_force5N_0829_20']
};

T_5N=[];
% Calculate the delta EMG for the monitorised channel 
for window = 2:length(tWindowInterval)
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T_5N = [T_5N; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
end
%%
%strapping force verification plotting
% a big table T_20%
T_20 = [T_2N(11:20,:);T_2N(31:40,:);T_5N(11:20,:);T_5N(31:40,:)];

F1 = figure(1);
boxchart(T_20.EMGChan, T_20.gradient,'GroupByColor',T_20.Amplitude);
ylabel('gradient ∆EMG'' [%]')
xlabel('Strapping Force - 5N (LHS) VS 2N (RHS)')
%ylim([-10,10]);
grid
set(gca,'FontSize',24)
F1.Position=[172 758 902 421];
F1.Children.XAxis.TickValues=[1:1:12];
legend('medium','high');
%F1.Children(2).Children(4).BoxFaceColor=[0,0,0];
%F1.Children(2).Children(3).BoxFaceColor=[0.00,0.34,0.74];
F1.Children(2).Children(2).BoxFaceColor=[0.93,0.69,0.13];
F1.Children(2).Children(1).BoxFaceColor=[0.49,0.18,0.56];

%% two-stimulator model contact area updated pilot test
nreps  = 5;                     % Number of repetitions of the trial
nChans = 14;                    % Number of channels EMG
Amp = [0, 85, 170, 255];        % Amplitudes considered
Optfreq = 120;                  % optimal frequency
tWindowInterval=[               % Window intervals to calculate delta EMG
    4,5;
    
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
    14,15;
    
    %{
    15,16;
    16,17;
    17,18 %after stimulation
    %}
    ];

Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/1207_xinyao_20_one']
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
    
    nchan = 2;
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_1 = [T0_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
    
end

disp(T0_1)
Path={
['/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/DataSample/1207_xinyao_20_two']
};

T0_2=[];
% Calculate the delta EMG for each channel and window
for window = 2:length(tWindowInterval)
    %{
    for nchan = 1:nChans
        fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
        T0_1 = [T0_1; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
    end
    %}
    
    nchan = 2;
    fprintf(strcat('Analyzing window interval: ',num2str(tWindowInterval(window,1)),'<t<',num2str(tWindowInterval(window,2)),', channel:',num2str(nchan),'...\n'))
    T0_2 = [T0_2; deltaEMG_2_2(nreps,[tWindowInterval(1,:);tWindowInterval(window,:)],Amp,0,nchan,Path{1})];
    
end

disp(T0_2)
%% Plotting pilot test (above)

window = '14<t<15';
%{
yData       = T.deltaEMG(strcmp(T.Period,window));
xData       = T.EMGChan(strcmp(T.Period,window));g
ColorGroups = T.Amplitude(strcmp(T.Period,window));
%}
yData       = T0_1.gradient;
xData       = T0_1.EMGChan;
ColorGroups = T0_1.Amplitude;
yData_2       = T0_2.gradient;
xData_2    = T0_2.EMGChan;
ColorGroups_2 = T0_2.Amplitude;

F1 = figure(1);
subplot(1,2,1);
boxchart(xData, yData,'GroupByColor',ColorGroups);
ylabel('gradient ∆EMG'' [%]')
xlabel('One-Stimulator Model')
ylim([-4,4]);
grid
set(gca,'xtick', [], 'FontSize',20)
subplot(1,2,2);
boxchart(xData_2, yData_2,'GroupByColor',ColorGroups_2);
ylabel('gradient ∆EMG'' [%]')
xlabel('Two-Stimulator Model')
ylim([-4,4]);
grid
legend('zero','low','medium','high');
set(gca,'xtick', [], 'FontSize',20)
F1.Position=[172 758 902 421];
F1.Children.XAxis.TickValues=[1:1:12];
F1.Children(2).Children(4).BoxFaceColor=[0,0,0];
F1.Children(2).Children(3).BoxFaceColor=[0.00,0.34,0.74];
F1.Children(2).Children(2).BoxFaceColor=[0.93,0.69,0.13];
F1.Children(2).Children(1).BoxFaceColor=[0.49,0.18,0.56];

%% Gradient over the channels

window = '14<t<15';
%{
yData       = T.deltaEMG(strcmp(T.Period,window));
xData       = T.EMGChan(strcmp(T.Period,window));
ColorGroups = T.Amplitude(strcmp(T.Period,window));
%}
yData       = T0_1.gradient;
xData       = T0_1.EMGChan;
ColorGroups = T0_1.Amplitude;

F1 = figure(1);
boxchart(xData, yData,'GroupByColor',ColorGroups);
ylabel('gradient ∆EMG'' [%]')
%ylim([-10,10]);
grid
set(gca,'FontSize',24)
F1.Position=[172 758 902 421];
F1.Children.XAxis.TickValues=[1:1:12];
legend('zero','low','medium','high');
%F1.Children(2).Children(4).BoxFaceColor=[0,0,0];
F1.Children(2).Children(3).BoxFaceColor=[0.00,0.34,0.74];
F1.Children(2).Children(2).BoxFaceColor=[0.93,0.69,0.13];
F1.Children(2).Children(1).BoxFaceColor=[0.49,0.18,0.56];
%% Delta EMG over the channels

window = '14<t<15';
%{
yData       = T.deltaEMG(strcmp(T.Period,window));
xData       = T.EMGChan(strcmp(T.Period,window));
ColorGroups = T.Amplitude(strcmp(T.Period,window));
%}
yData       = T0_1.deltaEMG;
xData       = T0_1.EMGChan;
ColorGroups = T0_1.Amplitude;

F1 = figure(1);
boxchart(xData, yData,'GroupByColor',ColorGroups);
ylabel('norm. ∆EMGrms [%]')
ylim([-10,10]);
grid
set(gca,'FontSize',24)
F1.Position=[172 758 902 421];
F1.Children.XAxis.TickValues=[1:1:12];
legend('zero','low','medium','high');
F1.Children(2).Children(4).BoxFaceColor=[0,0,0];
F1.Children(2).Children(3).BoxFaceColor=[0.00,0.34,0.74];
F1.Children(2).Children(2).BoxFaceColor=[0.93,0.69,0.13];
F1.Children(2).Children(1).BoxFaceColor=[0.49,0.18,0.56];

%% Raincloud per channel

% Channel to plot
chan  = 4;
% Fig. properties
fsize=31;


cb=[                    % Predefine colors for the plot
    0,0,0;
    0,0.45,0.74;
    0.93,0.69,0.13;
    0.49,0.18,0.56;
    0.47,0.67,0.19;
    0.3,0.75,0.93;
    0.64,0.08,0.18
    ];
f1=figure;
for ii = 1 : length(Amp)
    % Data to plot
    yData = T0_1.deltaEMG( T0_1.Amplitude == Amp(ii) & T0_1.EMGChan == chan );
    
    % Plotting
    h{ii} = raincloud_plot(yData',...
        'box_on', 1,...
        'box_dodge', 1,...
        'box_dodge_amount',1.5*ii,...
        'dot_dodge_amount', 1.5*ii,...
        'color', cb(ii,:),...
        'alpha',0.5,...
        'cloud_edge_col', cb(ii,:));
    box off;
end

% Figure properties
xlimit = [-10,10];
set(gca, 'XLim', xlimit);
set(gca, 'YTickLabel',[],'YTick',[]);
set(gca,'Color','none');
set(gca,'FontSize', fsize);     
legend([h{1}{1} h{2}{1} h{3}{1} h{4}{1}], {'zero','low','medium','high'});
f1.Position=[27 534 1031 485];
f1.Position=[27 497 748 522];
f1.Children(2).YLimMode = 'auto';
f1.Children(2).XTickLabelRotation = 270;
f1.Children(2).XAxisLocation = 'top';
f1.Children(2).XGrid = 'on';
f1.Children(2).Title.String = 'norm. ∆EMGrms [%]';
f1.Children(2).Title.FontWeight = 'normal';

%% Polar plot rms of mean all trials. - Fast way to observe the results

figure(2);
color(1,:) = [0,0,0];
color(2,:) = [0.00,0.34,0.74];
color(3,:) = [0.93,0.69,0.13];
color(4,:) = [0.58,0.76,0.93];

PolarPlot=figure(2);
for j = 1:length(Amp)
    rmsAEMG = zeros(1,nChans);
    stdAEMG = zeros(1,nChans);
    for i = 1:nChans
        rmsAEMG(i) = nanmean(T0_1.deltaEMG(T0_1.EMGChan==i & T0_1.Amplitude==Amp(j)));
        stdAEMG(i) = nanstd( T0_1.deltaEMG(T0_1.EMGChan==i & T0_1.Amplitude==Amp(j)));
    end
    theta = deg2rad([(0):30:(330)]);
    PolarPlot = polarplot(theta,rmsAEMG,'-o','linewidth',2,'color',color(j,:));
    hold on
end

% Plot properties
for i=1:12
     PolarPlot.Parent.ThetaTickLabel{i}=i;
end
PolarPlot.Parent.RLim = [-10,10];
PolarPlot.Parent.Parent.CurrentAxes.RAxis.Label.String = 'norm. ∆EMGrms [%]';
legend('zero','low','medium','high')
set(gca,'FontSize',20)

