%
%processes raw signal, extracts spikes. loads data directly from neuralynx
%file.
%
%filename: full path/name to neuralynx file
%totNrSamples: nr of samples in this file
%Hd: filter specification
%howManyBlocks: >0 do only so many blocks
%               0 : as many as there are
%startWithBlock: 1=start with first block
%                >1: dont process first x blocks
%includeRange: from/to timestamps of periods to sort. if empty everything
%is taken.
%
%extractionThreshold: how many times the running std of the energy signal
%to extract ?
%
%prewhiten: yes/no
%
%returns:
%blockOffsets -> start of each processed block in absolute timestamps
%
%
%urut/april04
%updated urut/feb07
%

% HM Edit
% Added spikeTimestamps as additional output - indices of timestamp for
% each spike peak

function [allSpikes, allSpikesNoiseFree, allSpikesCorrFree, allSpikesTimestamps, dataSamplesRaw,filteredSignal, rawMean,rawTraceSpikes,runStd2,upperlim,stdEstimates, blocksProcessed, noiseTraces, dataSamplesRawUncorrected, blockOffsets, allSpikeInds] = processRaw(filename, totNrSamples, Hd, params )
blockOffsets=[];

howManyBlocks = params.howManyBlocks;
startWithBlock = params.startWithBlock;
includeRange = params.includeRange;
%extractionThreshold = params.extractionThreshold;
%prewhiten = params.prewhiten;
alignMethod = params.alignMethod;

params.nrNoiseTraces=50;

%detect spikes
[runs, blocksize] = determineBlocks( totNrSamples,params.blocksize);

allSpikes=[];
allSpikesTimestamps=[];
allSpikesNoiseFree=[];
allSpikesCorrFree=[];
allSpikeInds = [];
allShifted = [];

%returns data of last block processed -- for debugging purposes
dataSamplesRaw=[];
dataSamplesRawUncorrected=[];
rawMean=[];
filteredSignal=[];
rawTraceSpikes=[];
runStd2=[];
upperlim=[];
stdEstimates=[];

noiseTraces=[];

if runs<startWithBlock
    disp(['warning: startWithBlock<runs']);
end

for i=startWithBlock:runs
    fromInd=(i-1)*blocksize+1;
    tillInd=i*blocksize;
    if tillInd>totNrSamples
        tillInd=totNrSamples;
    end

    if params.rawFileVersion<=2 || params.rawFileVersion==6
        disp(['CSC run: ' num2str(i) ' of ' num2str(runs) ' ind from ' num2str(ceil(fromInd/512)) ' to ' num2str(tillInd/512)]);
    else
        disp(['TXT run: ' num2str(i) ' of ' num2str(runs) ' ind from ' num2str(fromInd) ' to ' num2str(tillInd) ]);
    end
        
    %load data from file
    t1=clock;
    [timestampsRaw,dataSamplesRaw] = getRawData( filename, fromInd, tillInd, params.rawFileVersion, params.samplingFreq );
    
    disp(['time for raw read ' num2str(etime(clock,t1))]);

    if params.doGroundNormalization==1
        disp( ['normalizing ground to channels ' num2str(params.normalizationChannels) ] );
        
        rawFilePrefix = copyFieldIfExists( params, 'prefix', 'A' );   %file prefix default is Axx.ncs if not set
        
        [ commonGround ] = computeNormalizedGround(params.pathRaw, params.normalizationChannels, fromInd, tillInd, rawFilePrefix, params.rawFilePostfix);
        dataSamplesRawUncorrected = dataSamplesRaw;	
        dataSamplesRaw = dataSamplesRaw - commonGround;
    else
        disp('dont normalize ground');
    end
    
    wasPartialBlock=0; %flag to indicate that block was included that was only partially part of included time range
    if size(includeRange,1)>0
        includeMask=zeros(1, length(timestampsRaw));
        for j=1:size(includeRange,1)
            inds = find( timestampsRaw >= includeRange(j,1) & timestampsRaw <= includeRange(j,2) );
            if length(inds)>0
                includeMask(inds) = 1;
            end
        end
        
        if sum(includeMask)==0
            disp(['all samples in this block are excluded,skip from:' num2str(timestampsRaw(1),10) ' to:' num2str(timestampsRaw(end),10)]);
            continue;
        else
            if sum(includeMask)<length(timestampsRaw)
                %only some are excluded
                disp(['some samples are excluded,but still include everything. from:' num2str(timestampsRaw(1),10) ' to:' num2str(timestampsRaw(end),10)]);
            
                wasPartialBlock=1;
            else
                disp(['all samples included. from:' num2str(timestampsRaw(1),10) ' to:' num2str(timestampsRaw(end),10)]);
            end
        end
    end

    t1=clock; 
    [rawMean, filteredSignal, rawTraceSpikes,spikeWaveforms, spikeTimestamps, runStd2, upperlim, noiseTracesTmp] = extractSpikes( dataSamplesRaw, Hd, params );
    disp(['time for extraction ' num2str(etime(clock,t1))]);

    stdEstimates(i) = std(filteredSignal);

    %save(['c:\temp\traces\traces_B' num2str(i) '.mat'],'filteredSignal','spikeWaveforms','noiseTracesTmp');
    
    %add to global storage of noise traces, adding block # into first
    %column
    if size(noiseTracesTmp,1)>1
        noiseTracesTmp2=[];
        noiseTracesTmp2(1:size(noiseTracesTmp,1) ,1) = ones( size(noiseTracesTmp,1), 1 )*i;
        noiseTracesTmp2(1:size(noiseTracesTmp,1),2:size(noiseTracesTmp,2)+1)=noiseTracesTmp;
        noiseTraces = [noiseTraces; noiseTracesTmp2];
    end
    
    %upsample and classification
    disp(['upsample and classify run ' num2str(i) ' nr spikes ' num2str(size(spikeWaveforms,1))]);

    
    %classify
    %[spikeWaveformsNegative, spikeTimestampsNegativeTmp,nrSpikesKilledNegative,killedNegative] = classifySpikes(-1 * spikeWaveformsNegative,spikeTimestampsNegative);
    
    %only denoise if enough traces available
    
    nrNoiseTraces=size(noiseTraces,1);
    disp(['nrNoiseTraces ' num2str(nrNoiseTraces)]);

    %store orig waveforms for later whitening. 
	spikeWaveformsOrig = spikeWaveforms; %before upsampling
        
    %upsample and re-align
    spikeWaveforms=upsampleSpikes(spikeWaveforms);
%     spikeWaveforms = realigneSpikes(spikeWaveforms, spikeTimestamps, alignMethod, stdEstimates(i),filteredSignal);  %3==type is negative, do not remove any if too much shift
    [spikeWaveforms newTimeStamps shifted indsToKeep] = realigneSpikes(spikeWaveforms, spikeTimestamps, alignMethod, stdEstimates(i),filteredSignal);  %3==type is negative, do not remove any if too much shift
    
    % HM Edit
    % Make sure only usable spikes are kept
    spikeWaveformsOrig = spikeWaveformsOrig(indsToKeep',:);
    spikeTimestamps = newTimeStamps;

    %convert timestamps
    spikeTimestampsTmp = spikeTimestamps;
    if length( spikeTimestampsTmp ) > 0
        spikeTimestampsConverted = convertTimestamps( timestampsRaw, spikeTimestampsTmp, params.samplingFreq, params.rawFileVersion );
    else
        spikeTimestampsConverted = [];
    end
    blockOffsets = [blockOffsets timestampsRaw(1)];
    
    %-- exclude spikes outside of included time range (from partial blocks
    %with overlap)
    if wasPartialBlock
    
        includeMaskSpikes  = zeros(1, length(spikeTimestampsConverted));
        
        
        for j=1:size(includeRange,1)
            inds = find( spikeTimestampsConverted >= includeRange(j,1) & spikeTimestampsConverted <= includeRange(j,2) );
            includeMaskSpikes(inds)=1;
        end
        
        % only keep timestamps and waveforms of spikes that were included
        indsUse = find(includeMaskSpikes);
        spikeTimestampsConverted = spikeTimestampsConverted(indsUse);
        spikeWaveforms = spikeWaveforms(indsUse,:);
        spikeWaveformsOrig = spikeWaveformsOrig(indsUse,:);
        spikeTimestamps = spikeTimestamps(indsUse);
        
        disp(['Exclude spikes not in includeRange: ' num2str(length(indsUse)) ' of ' num2str(length(includeMaskSpikes)) ' included.']);
    end
    
    %put into global data structure
    allSpikes = [allSpikes; spikeWaveforms];
    allSpikesTimestamps = [allSpikesTimestamps spikeTimestampsConverted];
    allSpikesCorrFree = [ allSpikesCorrFree; spikeWaveformsOrig]; 
    allSpikeInds = [allSpikeInds (fromInd+spikeTimestamps-1)];
    allShifted = [allShifted shifted];

    if length(allSpikesTimestamps) ~= size(allSpikes,1)
        warning('nr timestamps is not equal to nr waveforms');
    end
    
    if length(allSpikeInds) ~= size(allSpikes,1)
        warning('nr timestamp indices is not equal to nr waveforms');
    end
    
    disp(['== number of spikes found in run ' num2str(i) ' ' num2str( size(spikeWaveforms,1)) ' total ' num2str(size(allSpikes,1)) ]);
    
    %stop as soon as 1000 spikes were detected,positive and negative
    if howManyBlocks==9999
        if size(allSpikes,1)>=1000 
            %&& size(allSpikesNegative,1)>=1000
	    warning('early stop is enabled - stopping detection because 1000 spikes detected.');
            blocksProcessed=i;
            break;
        end
    end
    
    %stop early if limit is set (by caller of function)
    if howManyBlocks>0
        if i>=howManyBlocks            
            blocksProcessed=i;
	    disp('early stop is enabled - stop detection because of tillBlocks limit');
            break;
        end
    end
end

blocksProcessed=i;

%whiten
if size(noiseTraces,1)>0 & size(allSpikesCorrFree,1)>0
    [trans, transUp, corr, stdWhitened] = posthocWhiten(noiseTraces, allSpikesCorrFree, alignMethod); % HM Edit - edited posthocWhiten to only upsample, without realigning
    allSpikesCorrFree = transUp;   
    
    % HM Edit
    if size(allSpikes,1) ~= size(allSpikesCorrFree,1)
        error('Spike waveforms and Spike waveform (corr free) don''t match in number');
    end
    for ii = 1:size(allSpikesCorrFree,1)
        currentSpike = allSpikesCorrFree(ii,:);
        diff = abs(allShifted(ii));
        maxPos = 95;
        if allShifted(ii)<0
            currentSpike = [currentSpike(1)*ones(1,diff) currentSpike(1:maxPos) currentSpike(maxPos+1:end-diff)]; % HM Edit
        elseif allShifted(ii)>0
            currentSpike = [currentSpike(diff+1:maxPos) currentSpike(maxPos+1:end) currentSpike(end)*ones(1,diff)];
        end
        allSpikesCorrFree(ii,:) = currentSpike;
    end
    
    %only store the autocorrelation,not all noise traces
    noiseTraces=corr;
end
