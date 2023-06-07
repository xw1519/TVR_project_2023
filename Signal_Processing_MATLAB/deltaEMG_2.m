% BLOCK I : ANALYSIS OPTIMAL STIMULATION AMPLITUDE
% A_OptimalFreq(ntrials,tWindowInterval,Amps,isPlot,Chan,path) 
% Analyse data from the BLOCK-1 for optimal stimulation amplitude
%
% Varibles:
% ntrials           - number of trials for each amplitude
% tWindowInterval   - windows to compare 
% Amps              - Stimulation amplitudes
% isPlot            - Plot the different trials
% Chan              - EMG channel to analyse
% path              - Path

    
function T=deltaEMG_2(ntrials,tWindowInterval,Amps,isPlot,Chan,path)

% Default values
    switch nargin
        case 0
            ntrials=5;
            tWindowInterval=[3,5;10,15];
            Amps=[0, 85, 170, 255];
            isPlot=1;
            Chan=3;
            path=uigetdir();

        case 1
            tWindowInterval=[3,5;10,15];
            Amps=[0, 85, 170, 255];
            isPlot=1;
            Chan=3;
            path=uigetdir();

        case 2
            Amps=[0, 85, 170, 255];
            isPlot=1;
            Chan=3;
            path=uigetdir();

        case 3
            isPlot=1;
            Chan=3;
            path=uigetdir();

        case 4
            Chan=3;
            path=uigetdir();

        case 5
            path=uigetdir();            
    end

    
    %% Load the path and define initial variables
    % Get all data files
    ToScan = strcat(path,'/*_QUAT*.mat');%'in windows: \'
    files = dir(ToScan);
    SortedFiles = natsortfiles({files.name});
    % Load random file to detect the number of channels
    load(fullfile(path,SortedFiles{1}))
    nChannels = size(RecInfo.Data.rmsEMG,2);
    
    % TABLE VARIABLES
    % The output of the function is a table with the output values of the
    % analysis:
    
    ControlConditions = length(Amps);                % Conditions to test
    AEMG = zeros(ntrials,ControlConditions);          % Array to store the delta EMG
    WindowBlock = {};                                 % Array to store the labels of the analysed window
   
    % 0 85 170 255 representing baseline, low, medium, high amplitude
    % trial 1
    % trial 2
    Amp = zeros(ntrials,ControlConditions);          % Array to store the labels per trial and amplitude
    
    % AUX VARIABLES
    cF = ones(1,length(Amps));                       % Counter per analysed trial
    
    % PLOTTING VARIABLES
    f = {};                                           % Variable to organise the plottings
    %% Color for plotting
    color(1,:) = linspace(1,0,ntrials);
    color(2,:) = zeros(1,ntrials);
    color(3,:) = linspace(0,1,ntrials);
    
    %% Extract data from the files
    for nfile=1:length(SortedFiles)
        % Load the raw EMG from the BIN file
        BinName = strcat(SortedFiles{nfile}(1:end-3),'bin');               
        rawEMG = readBin_simple([nChannels,Inf],'int16',fullfile(path,BinName));
        % Load additional data from the *.mat file
        FileToLoad = fullfile(path,SortedFiles{nfile});
        load(FileToLoad);
        
        % Get trial stimulation amplitude
        VibAmp = RecInfo.Experiment.Order(nfile);
        
        % If manual procecssing from raw data
        %{
        rmsWindow = 500;
        BufSize = 40;
        windowCount = 1;
        rmsEMG = zeros(nChannels, round((length(rawEMG)/BufSize) - (rmsWindow/BufSize))) ;
        samplect = 1;
        
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
        aux = NrmsEMG; % If offline analysis of the EMG
        %}

        aux = RecInfo.Data.rmsNormEMG'; % Using the online EMG samples during the trial
        Data = aux(Chan,:);
        
        % Remove the time offset
        Time = RecInfo.Data.TimeStamp - ones(length(RecInfo.Data.TimeStamp),1) .* RecInfo.Data.TimeStamp(1);
        
        % Get start part to analyse
        SamplesToAnalyse = find(Time > tWindowInterval(1,1) & Time < tWindowInterval(1,2));
        StartEMG = rms( Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) );
            
        % Get end part to analyse
        SamplesToAnalyse = find(Time > tWindowInterval(2,1) & Time < tWindowInterval(2,2));
        EndEMG = rms( Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) );
       
        % Calculate the delta EMG -  difference between beggining and end
        % of the trial
        difEMG = EndEMG - StartEMG;

        % Find the index to store the results 
        Fs = find(Amps == VibAmp); 
        % Store the results on the defined position
        AEMG(cF(Fs),Fs) = difEMG;
        Amp(cF(Fs),Fs) = VibAmp;
        WindowBlock(cF(Fs),Fs) = {strcat(num2str(tWindowInterval(2,1)),'<t<',num2str(tWindowInterval(2,2)))};

        %Plotting
        if isPlot
            f{Fs}=figure(Fs);
            hold on
            minLength  = min (length(Time), length(Data));
            plot(Time(1:minLength),Data(1:minLength).*100,'color',color(:,cF(Fs)),'linewidth',2)
            ylabel('Norm. rms EMG [%]');
            xlabel('Time [s]');
            xlim([0,20]);
            ylim([0,30])
            set(gca,'fontsize',30);
            f{Fs}.Position=[175 630 1202 389];
        end
        
        % update counter per analysed trial and condition 
        cF(Fs) = cF(Fs) + 1;
    end

    %% Data arrangement
    % reformat the data to save it in a table
    deltaEMG    = reshape(AEMG,1,[])'.*100; 
    Amplitude   = reshape(Amp,1,[])';
    Period      = reshape(WindowBlock,1,[])';
    SubjectID   = reshape(repmat({RecInfo.Session.SubjectID},1,ControlConditions*ntrials),1,[])';
    EMGChan     = repmat(Chan,1,ControlConditions*ntrials)';

    % Gather the data in a table
    T=table(Period,Amplitude,deltaEMG,SubjectID,EMGChan);
    
    % Plot results
    CatAmp = categorical(Amplitude,[0, 85, 170, 255]);
    if isPlot
        f{length(Amps)+1} = figure(length(Amps)+1);
        boxchart(CatAmp,T.deltaEMG);
        ylabel('norm. ∆EMGrms [%]')
        ylim([0,10]);
        grid
        set(gca,'FontSize',24)
        f{length(Amps)+1}.Position = [745 758 329 421];
    end

   %{
    Plot 3D diagram of Amplitude as X axis + Frequency as Y axis 
    if isPlot
        f{length(Amps)+1} = figure(length(Amps)+1);%?? -> need 0-255 calibration 
        scatter3(CatAmp, Freq, T.deltaEMG, 30, c, 'filled')%Freq real-time detection?
        xlabel('Amplitude')
        xlim([0,1]);
        ylabel('Frequency [Hz]')
        ylim([130,150]);
        zlabel('norm. ∆EMGrms [%]')
        zlim([0,10]);
        grid
        set(gca,'FontSize',24)
        f{length(Amps)+1}.Position = [745 758 329 421];
     end 
    %}

end
