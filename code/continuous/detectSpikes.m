%Spike detection and extraction; extraction of noise traces
%
%detects spikes from raw signal; re-aligns spikes. supports multiple methods of detection as well as alignment.
%this function can detect spikes in two ways: either by thresholding an arbitrary signal (thresholdSignal) or by picking spikes
%at pre-defined timepoints (searchInds). Also, the function supports different ways of finding the peak location. Some of them
%require a separate signal for this (peakFindSignal).
%
%inputs
%======
%rawSignal: bandpass filtered signal
%realRawMean: running mean of unfiltered signal
%thresholdSignal:signal to threshold for detection. empty if detection method used doesnt require thresholding.
%runningThres:threshold to use on runningStd. empty if detection method used doesnt require thresholding.
%params:struct,with fields: 
%       samplingFreq: sampling freq of rawSignal
%       detectionMethod: which method is used for detection (see extractSpikes.m)
%       limit: absolute value of maximal valid real signal. higher/lower then this is considered out of band.
%       peakAlignMethod: which method is used for peak detection (see extractSpikes.m)
%       alignMethod: additional parameters for peak finding. only used if peakAlignMethod==1 (see findPeak.m for values)
%       nrNoiseTraces: 0 (no) or >0 : get maximally x noise traces,each the same length as a waveform
%searchInds: localization of peaks, empty if thresholding is used as detection method. (depends on params.detectionMethod)
%peakFindSignal: signal used for peak finding. empty if peak finding method doesnt require this.
%
%
%outputs
%=======
%rawTrace: raw trace with only the parts left which were taken as spikes
%spikeWaveforms: raw waveforms
%spikeTimestamps: raw timestamp values, at peak
%noiseTraces: noise traces of specified length
%
%orig: urut/2004
%
%updated: 
%urut/feb07: parameterized dynamic range/out of range limit.
%urut/april07: allow direct passing of search Inds as alternative to to-be-thresholded signal.
%urut/may07: implemented new peak alignment methods (power,MTEO).
%
%===========================
function [rawTrace, spikeWaveforms, spikeTimestamps,  noiseTraces] = detectSpikes(rawSignal, realRawMean, thresholdSignal, runningThres, params, searchInds, peakFindSignal )
defineSortingConstants;

if nargin<=5
    searchInds=[];
    powerSignal=[];
end

switch ( params.samplingFreq )
    case 24000
        rawTraceLength=64;
        beforePeak=24;
        afterPeak=39;
    case 25000
        rawTraceLength=64;
        beforePeak=24;
        afterPeak=39;
    case 32556
        rawTraceLength=84;        
        beforePeak=24;
        afterPeak=59;
    case 32552
        rawTraceLength=84;        
        beforePeak=24;
        afterPeak=59;
    case 32000
        rawTraceLength=84;        
        beforePeak=24;
        afterPeak=59;        
    case 30000
        % HM Edit - later on, the spikes get upsampled to 256 points, but
        % if peak is left at 25, there will be a lot of realigning to do
        % (from ~76 to 95), so I've shifted the peak here to 31
        rawTraceLength = 84;
        beforePeak = 30;
        afterPeak = 53;
%         rawTraceLength=84;        
%         beforePeak=24;
%         afterPeak=59;        
    otherwise
        warning(['unknown sampling freq - assuming default values. Fs is ' num2str(params.samplingFreq) ] );
        rawTraceLength=64;
        beforePeak=24;
        afterPeak=39;
end

noiseTracesLength=rawTraceLength; 
        
spikeWaveforms=[];
spikeTimestamps=[];
noiseTraces=[];

stdRawSignal = std(rawSignal);

if params.detectionMethod<5   %length(searchInds)==0
    %find where signal crosses threshold (on power)
    foundind = ( thresholdSignal > runningThres );

    %convolute to make broader on both sides
    S1 =  filter(  ones(1,40), 1 , foundind)>=1;
    S2 =  filter(  ones(1,40), 1 , fliplr(foundind))>=1;

    inds=S1+S2;

    %--dont take parts where signal was out of band
    %limit=2046;
    take = realRawMean > -1*params.limit & realRawMean < params.limit;

    inds = inds & take;
    searchInds = find( inds >= 1);
end


maxCovered=0;    
counter=1;
counterPos=1;
totLength=length(rawSignal);

covered=[0 0];

rawTrace=zeros(1,length(rawSignal));

toInd2=0;

t1=clock;
notSig=0;

% HM commented out
% for i=1:length(searchInds)
%         if maxCovered >= searchInds(i)
%             continue;
%         end
%         
%         fromInd=searchInds(i)-beforePeak;
%         toInd=searchInds(i)+afterPeak;
%         if fromInd<=0 || toInd>totLength || fromInd<maxCovered
%             continue;
%         end
% 
%         spikeSignal = rawSignal( fromInd:toInd );
%         
%         switch ( params.peakAlignMethod )
%             case METHOD_PEAKFIND_FINDPEAK %determine peak with standard method
%                 peakInd=findPeak(spikeSignal, stdRawSignal, params.alignMethod);
%             case METHOD_PEAKFIND_NONE %no peak finding
%                 if params.detectionMethod==5  %wavelets determine their own peak time.
%                     peakInd=beforePeak;
%                 else
%                     peakInd=0;
%                 end
%             case METHOD_PEAKFIND_POWER %use power to find peak
%                 pSig = peakFindSignal(fromInd:toInd);                
%                 peakInd=findPeakPower(spikeSignal, pSig);
%             case METHOD_PEAKFIND_MTEO %use MTEO signal to find peak
%                 pSig = peakFindSignal(fromInd:toInd);                
%                 peakInd=findPeakMTEO(spikeSignal, pSig);
%             otherwise
%                 error('unknown peak detection method');
%         end
%         
%         %unclear spike, ignore it.
%         if peakInd==-1
%             notSig=notSig+1;
%             maxCovered=toInd;
%             %['not sig ' num2str(notSig)]
%             continue;
%         end
%         
% %         fromInd2 = fromInd+peakInd-beforePeak;
%         fromInd2 = fromInd+peakInd-beforePeak-1; % HM Edit
% %         toInd2   = fromInd+peakInd+afterPeak;
%         toInd2   = fromInd+peakInd+afterPeak-1; % HM Edit
%         %already covered,to prevent repeats in any case
%         if fromInd2<=0 || toInd2>length(rawSignal) || length(find(covered(:,1)==fromInd2))>0  && length(find(covered(:,2)==toInd2))>0
%             continue;
%         end
%         
%         covered(counterNeg,1:2) = [fromInd2 toInd2];
%         spikeWaveforms(counterNeg,:) = rawSignal( fromInd2:toInd2 )' ;
% %         spikeTimestamps(counterNeg) = fromInd+peakInd;
%         spikeTimestamps(counterNeg) = fromInd+peakInd-1; % HM Edit
% 
%         counterNeg=counterNeg+1; %how many spikes extracted so far
%         rawTrace(fromInd2:toInd2)=rawSignal( fromInd2:toInd2 );
%         maxCovered=toInd2;
% end

% HM Edit - to rank spikes within segments that exceed threshold, and
% extract windows based on rank
diffSearchInds = diff(searchInds);
searchIndBorders = find(diffSearchInds>1);

for i=1:length(searchIndBorders)
        
        % Grab the signal segment that crosses threshold
        if i == 1
            fromInd = searchInds(1);
            toInd = searchInds(searchIndBorders(i));
        elseif i == size(searchIndBorders,1)
            fromInd = searchInds(searchIndBorders(end));
            toInd = searchInds(size(searchInds,1));
        else
            fromInd = searchInds(searchIndBorders(i-1)+1);
            toInd = searchInds(searchIndBorders(i));
        end

        spikeSignal = rawSignal( fromInd:toInd );
        
        % Find the local peaks (positive and negative) and rank them (HM Edit) 
        [peaks,peak_locs] = findpeaks(abs(spikeSignal));
        peaks_ranked = flipud(sortrows([peaks peak_locs],1)); % [amp, ind]
        
        % HM Edit - Use the following if detecting spikes using peak rank
        % Extract exclusive windows of signal around each peak according to
        % rank (HM Edit)
        segment_inds = fromInd:toInd;
        peakInd = beforePeak + 1;
        inds_overlap = [];
        for jj = 1:size(peaks_ranked,1)
            % Check if this spike is overlapping with already stored spikes
            if ismember(jj,inds_overlap)
                continue;
            end
            % Align spike
            fromInd2 = segment_inds(peaks_ranked(jj,2))-beforePeak;
            toInd2 = segment_inds(peaks_ranked(jj,2))+afterPeak;
            spike_window = fromInd2:toInd2;
            % Skip spike if cannot capture full window within this data block
            if spike_window(1) < searchInds(1) || spike_window(end) > searchInds(end) 
                continue;
            end
            % Store extracted spikes
            if abs(rawSignal(segment_inds(peaks_ranked(jj,2)))) > runningThres
                covered(counter,1:2) = [fromInd2 toInd2];
                spikeWaveforms(counter,:) = rawSignal( fromInd2:toInd2 )' ; 
                spikeTimestamps(counter) = fromInd2+peakInd-1; % HM Edit

                counter=counter+1; %how many spikes extracted so far
                rawTrace(fromInd2:toInd2)=rawSignal( fromInd2:toInd2 );
                maxCovered=toInd2;
                % Adjust list of peaks so there are no overlapping windows
                % (any other peak in a 45 point time window (1.5ms) is excluded)    
%             inds_overlap = [inds_overlap; find(peaks_ranked(:,2) < peaks_ranked(jj,2)+afterPeak & peaks_ranked(:,2) > peaks_ranked(jj,2)-beforePeak)];
                inds_overlap = [inds_overlap; find(peaks_ranked(:,2) < peaks_ranked(jj,2)+22 & peaks_ranked(:,2) > peaks_ranked(jj,2)-22)];
            end
            
        end
        
% Sort spikes in ascending temporal order
[spikeTimestamps,indsort] = sort(spikeTimestamps);
spikeWaveforms = spikeWaveforms(indsort,:);




        
%         % HM Edit - Use the following if detecting spikes with original OSort code
%         if maxCovered >= searchInds(i)
%             continue;
%         end
%         
%         % Grab the signal segment that crosses threshold
%         if i == 1
%             fromInd = searchInds(1);
%             toInd = searchInds(searchIndBorders(i));
%         elseif i == size(searchIndBorders,1)
%             fromInd = searchInds(searchIndBorders(end));
%             toInd = searchInds(size(searchInds,1));
%         else
%             fromInd = searchInds(searchIndBorders(i-1)+1);
%             toInd = searchInds(searchIndBorders(i));
%         end
%         
%         fromInd=searchInds(i)-beforePeak;
%         toInd=searchInds(i)+afterPeak;
%         if fromInd<=0 || toInd>totLength || fromInd<maxCovered
%             continue;
%         end

%         spikeSignal = rawSignal( fromInd:toInd );
%         switch ( params.peakAlignMethod )
%             case METHOD_PEAKFIND_FINDPEAK %determine peak with standard method
%                 peakInd=findPeak(spikeSignal, stdRawSignal, params.alignMethod);
%             case METHOD_PEAKFIND_NONE %no peak finding
%                 if params.detectionMethod==5  %wavelets determine their own peak time.
%                     peakInd=beforePeak;
%                 else
%                     peakInd=0;
%                 end
%             case METHOD_PEAKFIND_POWER %use power to find peak
%                 pSig = peakFindSignal(fromInd:toInd);                
%                 peakInd=findPeakPower(spikeSignal, pSig);
%             case METHOD_PEAKFIND_MTEO %use MTEO signal to find peak
%                 pSig = peakFindSignal(fromInd:toInd);                
%                 peakInd=findPeakMTEO(spikeSignal, pSig);
%             otherwise
%                 error('unknown peak detection method');
%         end
%         
%         %unclear spike, ignore it.
%         if peakInd==-1
%             notSig=notSig+1;
%             maxCovered=toInd;
%             %['not sig ' num2str(notSig)]
%             continue;
%         end
%         
% %         fromInd2 = fromInd+peakInd-beforePeak;
%         fromInd2 = fromInd+peakInd-beforePeak-1; % HM Edit
% %         toInd2   = fromInd+peakInd+afterPeak;
%         toInd2   = fromInd+peakInd+afterPeak-1; % HM Edit
        %already covered,to prevent repeats in any case
%         if fromInd2<=0 || toInd2>length(rawSignal) || length(find(covered(:,1)==fromInd2))>0  && length(find(covered(:,2)==toInd2))>0
%             continue;
%         end
%         
%         covered(counterNeg,1:2) = [fromInd2 toInd2];
%         spikeWaveforms(counterNeg,:) = rawSignal( fromInd2:toInd2 )' ;
% %         spikeTimestamps(counterNeg) = fromInd+peakInd;
%         spikeTimestamps(counterNeg) = fromInd+peakInd-1; % HM Edit
% 
%         counterNeg=counterNeg+1; %how many spikes extracted so far
%         rawTrace(fromInd2:toInd2)=rawSignal( fromInd2:toInd2 );
%         maxCovered=toInd2;
end



%-- end extracting spikes
t2=clock;

%-- extract noise traces, if spikes were found (done using segments that
%don't contain spikes, beginning from first ind)
noiseTraces=[];
if params.nrNoiseTraces>0 & counter>1
    
    %if detection method offers no inds,make one
    if ~exist('inds')
        inds=zeros(1,length(rawSignal));
        for i=1:size(covered,1)
           inds( covered(i,1):covered(i,2) ) = 1; 
        end
    end
    
    searchInds = find( inds == 0);
    c=0;
    currentInd=1;
    if length(searchInds)>noiseTracesLength*params.nrNoiseTraces
        taken=[];
        totLengthRawSignal=length(rawSignal);
        
        while (1)
            if c>=params.nrNoiseTraces
                break;
            end
            
            indFrom=searchInds(currentInd);
            indTo=indFrom+noiseTracesLength-1;

            %if it overlaps into a spike,skip
            if sum( inds(indFrom:indTo) )>0
                currentInd=currentInd+noiseTracesLength;
                continue;
            end

            %if end is after end of the signal
            if indTo > totLengthRawSignal
                break;
            end

            c=c+1;
            taken(c,1:2)=[indFrom indTo];
            currentInd = currentInd+10+noiseTracesLength;

            if currentInd>length(searchInds)
                break;
            end
        end

        if size(taken,1)>0
            noiseTraces=zeros(size(taken,1), noiseTracesLength);
            for jj=1:size(taken,1)
                noiseTraces(jj,:)=rawSignal(taken(jj,1):taken(jj,2))';
            end
        end
        
    end
end
%---end extracting noise traces

