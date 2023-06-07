function T=deltaEMG_3(ntrials,tWindowInterval,Amps,isPlot,Chan,path)
%% Load the path and define initial variables
    % Get all data files
    ToScan = strcat(path,'/*_QUAT*.mat');%'in windows: \'
    files = dir(ToScan);
    SortedFiles = natsortfiles({files.name});
    
    ControlConditions = length(Amps);                % Conditions to test
    AEMG = zeros(ntrials,ControlConditions);         % Array to store the delta EMG
    WindowBlock = {};                                % Array to store the labels of the analysed window
    log = zeros(ntrials,ControlConditions);          % Array to store all gradients
   
    Amp = zeros(ntrials,ControlConditions);          % Array to store the labels per trial and amplitude
    nChannels = 14;                                  % Predefine
    % AUX VARIABLES
    cF = ones(1,length(Amps));                       % Counter per analysed trial
    
    % PLOTTING VARIABLES
    f = {};                                          % Variable to organise the plottings

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
        
        % raw data processed with notch filter
        % Detect noise peak frequency
        
        [pxx2,f2] = eH(rawEMG(Chan,:),1000,1,0);
        pM = max(abs(pxx2));    %magnitude
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
                %c=c+1;
            end
        end
        % filter the noise at 120Hz
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
        % filter the noise at 120Hz
        Qfactor = 10;
        fe = 83; % Assign peak frequency
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
        % filter the noise at 120Hz
        Qfactor = 10;
        fe = 191; % Assign peak frequency
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
        %{
        % filter the noise at peak frequency ~ 50Hz
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
                %c=c+1;
            end
        end
        %}
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
        
        %aux = RecInfo.Data.rmsNormEMG'; % Using the online EMG samples during the trial
        Data = aux(Chan,:);
        
        % Remove the time offset
        Time = RecInfo.Data.TimeStamp - ones(length(RecInfo.Data.TimeStamp),1) .* RecInfo.Data.TimeStamp(1);

        % Calculate gradient for time window 5s to 15s
        SamplesToAnalyse = find(Time > 5 & Time < 15);

        % Get start part to analyse
        SamplesToAnalyse = find(Time > tWindowInterval(1,1) & Time < tWindowInterval(1,2));
        StartEMG = rms( Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) );
        y_EMG = Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) ;
        X_time = Time(SamplesToAnalyse(1):SamplesToAnalyse(end));
        b = polyfit(X_time,y_EMG, 1);
        b_gradient = b(1); 
            
        % Get end part to analyse
        SamplesToAnalyse = find(Time > tWindowInterval(2,1) & Time < tWindowInterval(2,2));
        EndEMG = rms( Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) );
        %y_EMG = Data(SamplesToAnalyse(1):SamplesToAnalyse(end)) ;
        %X_time = Time(SamplesToAnalyse(1):SamplesToAnalyse(end)); 
        %b_gradient_end = regress(y_EMG',X_time);
       
        % Calculate the delta EMG -  difference between beggining and end
        % of the trial
        difEMG = EndEMG - StartEMG;

        % Calculate the delta gradient -  difference between two windows -
        % positive means increasing activity
        %difgradient = b_gradient_end_sum - b_gradient_start_sum;

        % Find the index to store the results 
        Fs = find(Amps == VibAmp); 
        % Store the results on the defined position
        AEMG(cF(Fs),Fs) = difEMG;
        Amp(cF(Fs),Fs) = VibAmp;
        log(cF(Fs),Fs) = b_gradient;
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
    %Period      = reshape(WindowBlock,1,[])';
    %SubjectID   = reshape(repmat({RecInfo.Session.SubjectID},1,ControlConditions*ntrials),1,[])';
    %EMGChan     = repmat(Chan,1,ControlConditions*ntrials)';
    %gradient    = reshape(log,1,[])'.*100;

    % Gather the data in a table
    T=table(Amplitude,deltaEMG);
    
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
end

