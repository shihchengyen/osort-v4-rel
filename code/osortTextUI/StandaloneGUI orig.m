%
%main functions for unsupervised sorting. this function is called by
%either the UI or by a script.
%
%urut/may07
function StandaloneGUI(pathsOrig, filesToProcess, thres, normalizationChannels, paramsIn)
% which Nlx2MatCSC.m

starttime=clock;

timeSortingStats=[]; %i time nrSpikesSorted
timeDetectionStats =[];%i time duration(in blocks of 512000 samples)

%set prefix of raw data files and result files
rawFilePrefix = copyFieldIfExists( paramsIn, 'rawFilePrefix', 'A' );
processedFilePrefix = copyFieldIfExists( paramsIn, 'processedFilePrefix', 'A' );

pathOutOrig=pathsOrig.pathOut;
for kkk=1:length(filesToProcess) % production
    % for kkk=1:length(filesToProcess)    % change to for for debugging
    paths=pathsOrig;
    
    i = filesToProcess(kkk);
    stri = num2str(i);
    
    currentThresInd=kkk;
    
    handles=[];
    handles.correctionFactorThreshold=0;  %minimal threshold, >0 makes threshold bigger
    handles.paramExtractionThreshold=thres(currentThresInd);
    
    handles = copyStructFields( paramsIn, handles, {'minNrSpikes','blockNrRawFig','doGroundNormalization','rawFileVersion','detectionMethod','detectionParams','peakAlignMethod','displayFigures',...
        'runningAverageLength','blocksize'});
    
    %handles.blockNrRawFig = paramsIn.blockNrRawFig;
    %handles.doGroundNormalization = paramsIn.doGroundNormalization;
    %handles.rawFileVersion=paramsIn.rawFileVersion;
    %handles.detectionMethod=paramsIn.detectionMethod;
    %handles.detectionParams=paramsIn.detectionParams;
    %handles.peakAlignMethod=paramsIn.peakAlignMethod;
    
    %define file format dependent properties
    [samplingFreq, limit, rawFilePostfix] = defineFileFormat(paramsIn.rawFileVersion, paramsIn.samplingFreq);
    
    handles.samplingFreq = samplingFreq; %sampling freq of raw data
    handles.limit = limit; %dynamic range
    
    handles.pathRaw = paths.pathRaw;
    
    %define include range
    handles.includeFilename=[paths.timestampspath '/timestampsInclude.txt'];
    includeRange=[];
    if exist(handles.includeFilename,'file')==2
        includeRange = dlmread(handles.includeFilename);
        
        ['include range is set from ' handles.includeFilename]
        handles.includeRange=includeRange; % kaminskij information for ploting spike rate
        
    else
        warning(['include range is not set! file not found: ' handles.includeFilename]);
        handles.includeRange=[]; % kaminskij information for ploting spike rate
    end
    
    
    %find the channels used for normalization of this electrode
    if paramsIn.doGroundNormalization
        if size(normalizationChannels,1)==1
            handles.normalizationChannels = normalizationChannels;
        else
            electrodeInd = mapChannelToElectrode( i );
            
            normalizeThisChannel = 0;
            if ~isempty( electrodeInd )
                % check if at least one normalization channel is given for this electrode
                if ~isempty(find(normalizationChannels(2,:)==electrodeInd))
                    normalizeThisChannel=1;
                end
            end
            
            if normalizeThisChannel
                handles.normalizationChannels = setdiff( normalizationChannels(1, find( normalizationChannels(2,:) == electrodeInd ) ), paramsIn.groundChannels);
            else
                disp(['normalization channel not defined for this channel, dont normalize - ' stri]);
                handles.normalizationChannels = [];
                handles.doGroundNormalization=0;
            end
        end
    else
        handles.normalizationChannels = [];
    end
    
    paths.pathOut = [ pathOutOrig '/' num2str(handles.paramExtractionThreshold)];
    
    if exist(paths.pathOut,'dir')==0
        ['creating directory ' paths.pathOut]
        mkdir(paths.pathOut);
    end
    
    handles.rawFilename=[paths.pathRaw '/' rawFilePrefix stri rawFilePostfix];
    if paramsIn.doDetection
        
        if exist(handles.rawFilename,'file')==0
            ['file does not exist, skip ' handles.rawFilename]
            continue;
        end
        
        handles.rawFilePrefix = rawFilePrefix;
        handles.rawFilePostfix = rawFilePostfix;
        
        starttimeDetection=clock;
        handles = StandaloneInit( handles , paramsIn.tillBlocks, paramsIn.prewhiten, paramsIn.alignMethod(kkk), includeRange );
        timeDetection = abs(etime(starttimeDetection,clock))
        
        %timeDetectionStats(size(timeDetectionStats,1)+1,:) = [ i timeDetection handles.blocksProcessedForInit];
        
        handles.filenamePrefix = [paths.pathOut processedFilePrefix stri];
        storeSortResultFiles( [], handles, 2 , 2 );%2==no figures, 2=noGUI
    end
    
    starttimeSorting=clock;
    
    %if exist('doFindThreshold','var')==1
    %    if doFindThreshold && exist(handles.rawFilename,'file')>0
    %        handles.prefix=processedFilePrefix;
    %        handles.from=num2str(i);
    %        handles.includeRange = includeRange;
    %        
    %        handles = copyStructFields( paramsIn, handles, {'tillBlocks','prewhiten','alignMethod(kkk)'});
    %        %handles.howManyBlocks = paramsIn.tillBlocks;
    %        %handles.prewhiten = paramsIn.prewhiten;
    %        %handles.alignMethod = paramsIn.alignMethod(kkk);
    %        
    %        thresholds=[3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 ];
    %        produceThresholdFig(handles, [paths.pathFigs num2str(thres(currentThresInd)) '/'], thresholds);
    %    end
    %end
    
    if paramsIn.doSorting || paramsIn.doFigures || ~paramsIn.noProjectionTest
        handles.basepath=paths.pathOut;
        handles.prefix=processedFilePrefix;
        handles.from=stri;
        [handles,fileExists] = loadSortResultFiles([],handles, 2);
        handles.filenamePrefix=[paths.pathOut handles.prefix stri];
        
        if fileExists==0
            ['File does not exist: ' handles.filenamePrefix];
            continue;
        end
        
        if paramsIn.doSorting
            if size(handles.newSpikesNegative,1)>0
                starttimeSorting=clock;
                
                [handles] = sortMain( [], handles, 2, paramsIn.thresholdMethod  ); %2=no GUI
                
                timeSorting=abs(etime(starttimeSorting,clock))
                %timeSortingStats(size(timeSortingStats,1)+1,:) = [i timeSorting length(handles.assignedClusterNegative)];
                
                storeSortResultFiles(  [], handles, 2, 2  );%2==no figures
            else
                'nothing to sort (0 spikes)'
            end
        end
        
        handles.label=[ paths.patientID ' ' handles.prefix handles.from ' Th:' num2str(thres(currentThresInd))];
        handles.label = strrep(handles.label,'_',' ');
        disp(['producing figure for ' handles.label]);
        
        %clusters and PCA figures
        if paramsIn.doFigures
            outputPathFigs=[paths.pathFigs num2str(thres(currentThresInd)) '/'];
            produceFigures(handles, outputPathFigs, paramsIn.outputFormat, paramsIn.thresholdMethod);
            produceOverviewFigure(handles, outputPathFigs, paramsIn.outputFormat, paramsIn.thresholdMethod)
        end
        
        %projection test
        if ~paramsIn.noProjectionTest
            produceProjectionFigures(handles,[paths.pathFigs num2str(thres(currentThresInd)) '/'], paramsIn.outputFormat, paramsIn.thresholdMethod);
        end
        
    end
    
    
    % make raw plots at the end, so that sorting result is included in raw plots as colored spikes (cluster identity)
    if paramsIn.doRawGraphs && exist(handles.rawFilename,'file')>0
        handles.prefix = processedFilePrefix;
        handles.rawFilePrefix = rawFilePrefix;
        handles.rawFilePostfix = rawFilePostfix;
        handles.from=stri;
        handles.basepath=paths.pathOut;
        
        [handles,fileExists] = loadSortResultFiles([],handles, 2);
        produceRawTraceFig(handles, [paths.pathFigs num2str(thres(currentThresInd)) '/'], paramsIn.outputFormat);
        
    end
    
    
   % plot zoomed-in raw traces to examine waveforms 
   if isfield(paramsIn, 'doshortRawGraphs') && exist(handles.rawFilename,'file')>0
        if paramsIn.doshortRawGraphs~=0
            handles.prefix = processedFilePrefix;
            handles.rawFilePrefix = rawFilePrefix;
            handles.rawFilePostfix = rawFilePostfix;
            handles.from=stri;
            handles.basepath=paths.pathOut;
            [handles,fileExists] = loadSortResultFiles([],handles, 2);
            produceRawTraceFigSHORT(handles, [paths.pathFigs num2str(thres(currentThresInd)) '/'], paramsIn.outputFormat, paramsIn.doshortRawGraphs);
            
            
        end
   end
    
   % plot cross correlations between units
    if isfield(paramsIn, 'produceXcorrFigures')
        if paramsIn.produceXcorrFigures~=0
            XcorrlengthinMS=40;
            binsize=1;
            produceXcorrFigures(handles,[paths.pathFigs num2str(thres(currentThresInd)) '/'], paramsIn.outputFormat, paramsIn.thresholdMethod,XcorrlengthinMS,binsize);
            
        end
    end
    
    
end

etime(clock,starttime)
