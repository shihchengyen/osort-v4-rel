%
% Standalone_textGUI_prod1.m
%
% production code, sort epilepsy data on server


%% which files to sort
paths=[];

paths.basePath=['V:\dataRawEpilepsy\P48CS_031017\NO']; %add getuserdir in brackets if files are in C: ex. [getuserdir'\Documents\Practice OSort\P39CS_080515'];
paths.pathOut=[paths.basePath '/sort/'];
paths.pathRaw=[paths.basePath '/raw/'];
paths.pathFigs=[paths.basePath '/figs/']; %When creating figs after  practice merges, use /figs_merged/ directory instead of just /figs/ so nothing gets overwritten

paths.timestampspath=paths.basePath;             % if a timestampsInclude.txt file is found in this directory, only the range(s) of timestamps specified will be processed. 

paths.patientID='P48CS'; % label in plots

filesToProcess = [1:64, 192:207];  %which channels to detect/sort
noiseChannels  = [  ]; 

% which channels to ignore

groundChannels=[ ]; %which channels are ground (they are ignored)

doGroundNormalization = 0;
normalizeOnly = [  ]; %which channels to use for normalization (common average subtraction)

if exist('groundChannels') && ~doGroundNormalization
  filesToProcess=setdiff(filesToProcess, groundChannels);
end

%default align is mixed, unless listed below as max or min.
filesAlignMax=filesToProcess;
filesAlignMin=[];  % all channels are min when filesToProcess is added instead of brackets

extractionThreshold=5.002;

%will add 0.01 or 0.02 automatically if all are max/min

%% global settings
paramsIn=[];

paramsIn.rawFilePrefix='CSC';        % some systems use CSC instead of A
paramsIn.processedFilePrefix='A';

paramsIn.rawFileVersion = 6; %1 is analog cheetah, 2 is digital cheetah (NCS), 3 is txt file, 6 is Atlas. See defineFileFormat.m
paramsIn.samplingFreq = 30000; %only used if rawFileVersion==3

%which tasks to execute
paramsIn.tillBlocks = 999;  %how many blocks to process (each ~20s). 0=no limit.
paramsIn.doDetection =1;
paramsIn.doSorting = 1;
paramsIn.doFigures = 1;
paramsIn.noProjectionTest = 1;
paramsIn.doRawGraphs = 1;
paramsIn.doGroundNormalization=doGroundNormalization;

paramsIn.displayFigures = 0 ;  %1 yes (keep open), 0 no (export and close immediately); for production use, use 0 

paramsIn.minNrSpikes=50; %min nr spikes assigned to a cluster for it to be valid
                                                                                                                         
%params
paramsIn.blockNrRawFig=[ 10 15 20 25 ];
paramsIn.outputFormat='png';
paramsIn.thresholdMethod=1; %1=approx, 2=exact
paramsIn.prewhiten=0; %0=no, 1=yes,whiten raw signal (dont)
paramsIn.defaultAlignMethod=3;  %only used if peak finding method is "findPeak". 1=max, 2=min, 3=mixed
paramsIn.peakAlignMethod=1; %1 find Peak, 2 none, 3 peak of power, 4 MTEO peak
                        
%for wavelet detection method
%paramsIn.detectionMethod=5; %1 power, 2 T pos, 3 T min, 4 T abs, 5 wavelet
%dp.scalesRange = [0.2 1.0]; %in ms
%dp.waveletName='bior1.5'; 
%paramsIn.detectionParams=dp;
%extractionThreshold=0.1; %for wavelet method

%for power detection method
paramsIn.detectionMethod=1; %1 power, 2 T pos, 3 T min, 3 T abs, 4 wavelet
dp.kernelSize=18; 
paramsIn.detectionParams=dp;
%extractionThreshold = 5.5;  % extraction threshold

%automatically adjust threshold to indicate min/max
if isequal(filesToProcess,filesAlignMax) & isempty(filesAlignMin)
	extractionThreshold=extractionThreshold+0.001;
end
if isequal(filesToProcess,filesAlignMin) & length(filesAlignMax)==0
	extractionThreshold=extractionThreshold+0.002;
end

thres         = [repmat(extractionThreshold, 1, length(filesToProcess))];

%% execute
[normalizationChannels,paramsIn,filesToProcess] = StandaloneGUI_prepare(noiseChannels,doGroundNormalization,paramsIn,filesToProcess,filesAlignMax, filesAlignMin, normalizeOnly, groundChannels);
StandaloneGUI(paths, filesToProcess, thres, normalizationChannels, paramsIn);
