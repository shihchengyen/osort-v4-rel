%% which files to sort
paths=[];
paths.basePath='c:/dataLocalCopy/DL_022309/';
paths.pathOut=[paths.basePath 'sort/'];
paths.pathRaw=[paths.basePath '/rawLFP/'];
paths.pathFigs=[paths.basePath '/figs/'];

paths.patientID='P23S4';

filesToProcess=[  2   ];  %which files to detect/sort
noiseChannels=[ ];%24 32 35 36 39 40  ];

groundChannels=[8 9 22 32 50];%s2

doGroundNormalization=0;
normalizeOnly=[]; %which channels to use for normalization

if exist('groundChannels') && ~doGroundNormalization
  filesToProcess=setdiff(filesToProcess, groundChannels);
end

filesAlignMax=[ 2:7  ];
filesAlignMin=[  ];


%% extraction threshold
%extractionThreshold=0.1; %for wavelet method

extractionThreshold=5;
thres         =[repmat(extractionThreshold, 1, length(filesToProcess))];

%% global settings
paramsIn=[];

paramsIn.rawFileVersion=2; %1 is analog cheetah, 2 is digital cheetah (NCS), 3 is txt file.  determines sampling freq&dynamic range.
paramsIn.samplingFreq=24000; %only used if rawFileVersion==3

%which tasks to execute
paramsIn.tillBlocks=999;  %how many blocks to process (each ~20s). 0=no limit.
paramsIn.doDetection=0;
paramsIn.doSorting=0;
paramsIn.doFigures=1;
paramsIn.noProjectionTest=0;
paramsIn.doRawGraphs=0;
paramsIn.doGroundNormalization=doGroundNormalization;

paramsIn.minNrSpikes=20; %min nr spikes assigned to a cluster for it to be valid
                                                                                                                         


%params
paramsIn.blockNrRawFig=[ 90 100 120  ];
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

%for power detection method
paramsIn.detectionMethod=1; %1 power, 2 T pos, 3 T min, 3 T abs, 4 wavelet
dp.kernelSize=18; 
paramsIn.detectionParams=dp;
        
%% pre-processing. below this lines are no parameters.

%% --mapping of ground to electrode
normalizationChannels=[]; %maps channel to electrode. All that are listed are used for normalization.
if exist('noiseChannels')
    filesToProcess=setdiff(filesToProcess, noiseChannels);
end

if doGroundNormalization
    electrodeAssignment = [ 1:64; repmat(1,1,8) repmat(2,1,8) repmat(3,1,8) repmat(4,1,8) repmat(5,1,8) repmat(6,1,8) repmat(7,1,8) repmat(8,1,8) ]; %which wire is on which electrode

    %map channels to El number
    for j=1:length(normalizeOnly)
        indOfEl = find( electrodeAssignment(1,:)==normalizeOnly(j) );
        normalizationChannels(1:2,j) = [normalizeOnly(j); electrodeAssignment(2, indOfEl)];
    end
    
    %for j=1:length(filesToProcess)
    %    if length(find(normalizeOnly==filesToProcess(j)))==1
    %        normalizationChannels(1:2,j) = [filesToProcess(j); electrodeAssignment(2, j)];
    %    end
    %end

    excludeChannels=[]; %setdiff(filesToProcess, filesToProcess);

    if exist('noiseChannels')
        excludeChannels=[excludeChannels noiseChannels];
    end

    %exclude noise channels from the grand average.
    excludeInds=[];
    for i=1:length(excludeChannels)
        excludeInds = [ excludeInds find( normalizationChannels(1,:) == excludeChannels(i) ) ];
    end
    includeInds = setdiff( 1:size(normalizationChannels,2), excludeInds);
    normalizationChannels = normalizationChannels(:, includeInds);

    paramsIn.groundChannels=groundChannels;
end

%align method can be changed for each channe;; alignMethod is only used if peakAlignMethod=1 (findPeak)
paramsIn.alignMethod = repmat(paramsIn.defaultAlignMethod, 1, length(filesToProcess));
for i=1:length(filesToProcess)
	if length(find( filesAlignMax == filesToProcess(i) ))==1
		paramsIn.alignMethod(i) = 1;
	end
	if length(find( filesAlignMin == filesToProcess(i) ))==1
		paramsIn.alignMethod(i) = 2;
	end
end

%% execute
StandaloneGUI(paths, filesToProcess, thres, normalizationChannels, paramsIn);
