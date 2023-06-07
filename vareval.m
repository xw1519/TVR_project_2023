% Load mat dat a\
path=uigetdir();

% Get all data files
ToScan = strcat(path,'/*_QUAT*.mat');%'in windows: \'
files = dir(ToScan);
SortedFiles = natsortfiles({files.name});
% Load random file to detect the number of channels
load(fullfile(path,SortedFiles{1}))
nChannels = size(RecInfo.Data.rmsEMG,2);