%
%post-hoc whitening & upsampling of waveforms
%
%used in cases when online estimate of covariance is unstable
%
%returns:
%trans: transformed spikes
%corr:autocorrelation function of the noise
%stdWhitened: std of the whitened waveforms. if this is substantion > 0 -> electrode moved.
%
%
%urut/nov05
function [trans, transUp, corr, stdWhitened] = posthocWhiten(noiseTraces, origWaveforms, alignMethod )
nrSamplesPerSpike=size(origWaveforms,2);

noiseTraces=noiseTraces(:,2:nrSamplesPerSpike+1);


n=size(noiseTraces,1);
if n>10000
    n=10000;
end

noiseTraces=noiseTraces(1:n,:);
noiseTraces=noiseTraces';
noiseTraces=noiseTraces(:);

%estimate autocorrelation
corr=xcorr(noiseTraces,nrSamplesPerSpike,'biased');
corr=corr(nrSamplesPerSpike+1:end-1);

%whiten
C1=toeplitz(corr);
C1inv=inv(C1);
R1=chol(C1inv);
trans = origWaveforms * R1';

%stdWhitened = mean(std(trans));

stdWhitened=0; %no min/max here (to enforce strict re-alignment regardless of significance for pre-whitened data)

%upsample and re-align
transUp = upsampleSpikes(trans);
% transUp = realigneSpikesWhiten(transUp, [], alignMethod, stdWhitened);  % HM Edit
% HM Edit - realign this way if finding spikes with ranked peaks
shifted = zeros(1,size(transUp,1));
for jj = 1:size(transUp,1)
    spike = transUp(jj,:);
    % Prevent mistaking of peak at boundaries of block
    spike_trim = spike;
    spike_trim(1:90) = 0; spike_trim(100:end) = 0;
    ind_peak = find(abs(spike_trim) == max(abs(spike_trim)));
    diff = ind_peak-95;
    if diff>0
        transUp(jj,:) = [spike(1+diff:end) repmat(spike(end),1,diff)];
    elseif diff<0
        transUp(jj,:) = [repmat(spike(1),1,abs(diff)),spike(1:end-abs(diff))];
    end
    shifted(jj) = diff;
end



