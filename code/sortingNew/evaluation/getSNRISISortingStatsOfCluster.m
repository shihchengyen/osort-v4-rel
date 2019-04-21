%
%gets, for all clusters of a particular channel, the %ISI < 3ms for each cluster as well as the SNR of the mean waveform.
%originally in evalSNRAllCells.m, moved to this function march07.
%
%
%urut/march07
function [stats,nrSorted, allWaveforms] = getSNRISISortingStatsOfCluster(channel, sessionID, useNegative, assignedNegative, newSpikesNegative, newTimestampsNegative, stdEstimateOrig)
stats=[];
nrSorted=0;

allWaveforms = []; %mean/sd/se of the raw waveforms
for k=1:length(useNegative)
    clNr = useNegative(k);
    inds = find( assignedNegative == clNr );
    
    mWaveform = mean ( newSpikesNegative(inds,:) );
    sWaveform = std ( newSpikesNegative(inds,:) );
    seWaveform = sWaveform./sqrt(length(inds));
    
    allWaveforms(k).m = mWaveform;
    allWaveforms(k).sd = sWaveform;
    allWaveforms(k).se = seWaveform;
       
    SNR = calcSNR( mWaveform, stdEstimateOrig);

    SNRPeak = abs(max(mWaveform))./stdEstimateOrig;
    
    ISI = diff(newTimestampsNegative(inds));
    ISI = ISI./1000; %to ms

    ISI_inSec = ISI./1000; %to sec
    
    ISIbelowAbs  = length(find( ISI < 3.0 ));
    ISIbelowPerc = ISIbelowAbs*100/length(ISI);
    
    CV = calcCV(ISI_inSec);    
    
    ignoreMode = 1;   
    CV2 = calcCV2(ISI_inSec, ignoreMode); 
    CV2tmp = calcCV2(ISI_inSec, 0); 

    [CV2 CV2tmp]
    
    BI = calcBurstIndex(ISI_inSec);
    mISIAll=mean(ISI_inSec);
    
    entry = [ sessionID channel clNr length(inds) stdEstimateOrig SNR max(mWaveform) ISIbelowPerc ISIbelowAbs CV CV2 BI mISIAll SNRPeak];
    entryNr = size(stats,1)+1;

    nrSorted=nrSorted+length(inds);

    stats(entryNr,:) = entry;
end
