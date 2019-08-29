%
% RunOSort.m
%
% Demo file on how to run OSort on Atlas data.
%

% HM Edit
function RunOSort(pathstr);

% get directory separator for this platform
filesepchar = filesep();
% convert pathstr to titlestr
pstr = strsplit(pathstr,filesepchar);
% we are expecting pathstr to be a channel directory, so the last 4 items
% in pstr are what we need
psend = size(pstr,2);
titlestr = [pstr{psend-4}(1) pstr{psend-3} pstr{psend-2}([1 end-1:end]) pstr{psend-1}([1 end-1:end]) pstr{psend}([1 end-2:end])];

%% which files to sort
paths=[];

% paths.basePath=['/Volumes/Hippocampus/Data/picasso/' Date '/' Session];
paths.basePath = '.';    
% paths.pathOut=[paths.basePath '/sort'];
% paths.pathRaw=[paths.basePath '/' Array];
paths.pathRaw = paths.basePath;
% if a timestampsInclude.txt file is found in this directory, only the 
% range(s) of timestamps specified will be processed. 
paths.timestampspath = paths.basePath;             

filesToProcessStr = pstr{psend}(end-2:end); % HM Edit
% filesToProcessStr = '1'; % HM Edit
%which channels to detect/sort
filesToProcess = str2double(filesToProcessStr);  
% which channels to ignore
noiseChannels  = [   ]; 

%which channels are ground (ignore)
groundChannels=[]; 

doGroundNormalization=0;
%which channels to use for normalization
normalizeOnly=[]; 

if exist('groundChannels') && ~doGroundNormalization
  filesToProcess=setdiff(filesToProcess, groundChannels);
end

%default align is mixed, unless listed below as max or min.
filesAlignMax=[ ];
filesAlignMin=[ ];


%% global settings
paramsIn=[];

% some systems use CSC instead of A
paramsIn.rawFilePrefix='channel';        
paramsIn.processedFilePrefix='Ch';

% see defineFileFormat.m for options
paramsIn.rawFileVersion = 5; 
% only used if rawFileVersion==3
paramsIn.samplingFreq = 30000; 

%which tasks to execute
%how many blocks to process (each ~20s). 0=no limit.
paramsIn.tillBlocks = 0;  
paramsIn.doDetection = 1;
paramsIn.doSorting = 1;
paramsIn.doFigures = 1;
paramsIn.noProjectionTest = 0;
paramsIn.doRawGraphs = 1;
% zoomed-in raw graphs
paramsIn.doshortRawGraphs = 1; 
paramsIn.doGroundNormalization=doGroundNormalization;

% 1 yes (keep open), 0 no (export and close immediately); for production 
% use, use 0 
paramsIn.displayFigures = 0 ;  

%min nr spikes assigned to a cluster for it to be valid
paramsIn.minNrSpikes=100; 
                                                                                                                         
%params
% for which blocks will a raw figure be made
paramsIn.blockNrRawFig=[ 1 ];   
paramsIn.outputFormat='png';
% 1=approx, 2=exact
paramsIn.thresholdMethod=2; 
% 0=no, 1=yes,whiten raw signal (dont)
paramsIn.prewhiten=0; 
% only used if peak finding method is "findPeak". 1=max, 2=min, 3=mixed
paramsIn.defaultAlignMethod=3;  
%1 find Peak, 2 none, 3 peak of power, 4 MTEO peak
paramsIn.peakAlignMethod=1; 
                        
% % for wavelet detection method
% 1 power, 2 T pos, 3 T min, 4 T abs, 5 wavelet
% paramsIn.detectionMethod=5; 
% in ms
% dp.scalesRange = [0.2 1.0]; 
% bior1.5 is default, db2
% dp.waveletName='bior1.5'; 
% paramsIn.detectionParams=dp;
% for wavelet method
% extractionThreshold=0.2; 

% for power detection method
% 1 power, 2 T pos, 3 T min, 4 T abs, 5 wavelet
paramsIn.detectionMethod=1; 
dp.kernelSize=25; 
paramsIn.detectionParams=dp;
% extraction threshold
extractionThreshold = 8;  

thres = [repmat(extractionThreshold, 1, length(filesToProcess))];

%% Filenames

if isfield(dp,'kernelSize')
    foldername = [paths.pathRaw filesepchar 'oSort' filesepchar ...
        'detect' num2str(paramsIn.detectionMethod) ...
        'Thresh' num2str(extractionThreshold) 'kern' num2str(dp.kernelSize)];
elseif isfield(dp,'waveletName')
    foldername = [paths.pathRaw filesepchar 'oSort' filesepchar ... 
        'detect' num2str(paramsIn.detectionMethod) dp.waveletName ...
        'Thresh' num2str(extractionThreshold)];
end

% if isfield(dp,'kernelSize')
%     foldername = [paths.pathRaw filesepchar 'oSort' filesepchar 'detect' ...
%     	num2str(paramsIn.detectionMethod) 'peakType' ...
%     	num2str(paramsIn.peakAlignMethod) 'Align' ...
%     	num2str(paramsIn.defaultAlignMethod) 'Thresh' ...
%     	num2str(extractionThreshold) 'kern' num2str(dp.kernelSize)];
%     dp.kernelSize
% elseif isfield(dp,'waveletName')
%     foldername = [paths.pathRaw filesepchar 'oSort' filesepchar 'detect' ...
%     	num2str(paramsIn.detectionMethod) dp.waveletName 'peakType' ...
%     	num2str(paramsIn.peakAlignMethod) 'Align' ...
%     	num2str(paramsIn.defaultAlignMethod) 'Thresh' ...
%     	num2str(extractionThreshold)];
% end

paths.pathOut=[foldername filesepchar 'sort'];
paths.pathFigs=[foldername filesepchar 'figs'];
indDate = strfind(paths.basePath,[filesepchar 'session']);
% label in plots
% paths.patientID=['Picasso' paths.basePath(indDate-8:indDate-1) 'Ch' filesToProcessStr]; 
paths.patientID = titlestr;

%% execute
[normalizationChannels,paramsIn] = StandaloneGUI_prepare(noiseChannels,...
	doGroundNormalization,paramsIn,filesToProcess,filesAlignMax, ...
	filesAlignMin, normalizeOnly, groundChannels);

StandaloneGUI(paths, filesToProcess, thres, normalizationChannels, paramsIn);
