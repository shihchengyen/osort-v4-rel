%
% Standalone_textGUI_demo_new.m
%
% Demo file on how to run OSort on Atlas data.
%

% HM Edit
function [] = Standalone_textGUI_demo_new(Date,Session,Array,Channel);


%% which files to sort
paths=[];

paths.basePath=['/Volumes/Hippocampus/Data/picasso/' Date '/' Session];    
% paths.pathOut=[paths.basePath '/sort'];
paths.pathRaw=[paths.basePath '/' Array];
% if a timestampsInclude.txt file is found in this directory, only the 
% range(s) of timestamps specified will be processed. 
paths.timestampspath=paths.basePath;             

filesToProcessStr = Channel; % HM Edit
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
paramsIn.processedFilePrefix='P';

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
paramsIn.blockNrRawFig=[ 100 ];   
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
dp.kernelSize=18; 
paramsIn.detectionParams=dp;
% extraction threshold
extractionThreshold = 5;  

thres = [repmat(extractionThreshold, 1, length(filesToProcess))];

%% Filenames

if isfield(dp,'kernelSize')
    foldername = [paths.pathRaw '/' paramsIn.rawFilePrefix ...
    	filesToProcessStr '/oSort' '/detect' num2str(paramsIn.detectionMethod) ...
    	'peakType' num2str(paramsIn.peakAlignMethod) 'Align' ...
    	num2str(paramsIn.defaultAlignMethod) 'Thresh' ...
    	num2str(extractionThreshold) 'kern' num2str(dp.kernelSize)];
    dp.kernelSize
elseif isfield(dp,'waveletName')
    foldername = [paths.pathRaw '/' paramsIn.rawFilePrefix filesToProcessStr ...
    	'/oSort' '/detect' num2str(paramsIn.detectionMethod) dp.waveletName ...
    	'peakType' num2str(paramsIn.peakAlignMethod) 'Align' ...
    	num2str(paramsIn.defaultAlignMethod) 'Thresh' num2str(extractionThreshold)];
end
paths.pathOut=[foldername '/sort'];
paths.pathFigs=[foldername '/figs'];
indDate = strfind(paths.basePath,'/session');
% label in plots
paths.patientID=['Picasso' paths.basePath(indDate-8:indDate-1) 'Ch' filesToProcessStr]; 

%% execute
[normalizationChannels,paramsIn] = StandaloneGUI_prepare(noiseChannels,...
	doGroundNormalization,paramsIn,filesToProcess,filesAlignMax, ...
	filesAlignMin, normalizeOnly, groundChannels);

StandaloneGUI(paths, filesToProcess, thres, normalizationChannels, paramsIn);
