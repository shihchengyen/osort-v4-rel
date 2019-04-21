%
% get normalized waveforms and metrics for them, such as trough-to-peak time
%
% normalizes all waveforms to peak/through and flips to same orientation
%
% waves is a matrix of normlized waveforms
% d is the peak-to-through metric (Mitchell et al 2007 Neuron). Also referred to as Valley-to-peak by some. in ms.
% metricPos is the through and peak pos that was used
% w is the half-width metric for the valley (in ms)
% vpRatio is the valley-to-peak amplitude ratio
%
%urut/MPI/oct11
function [waves,d, metricPos,w, vpRatio ] = calcWaveformMetrics( peakPos, Fs, allWaveforms )
metricPos=[];
normMode=2; %0 none, 1 all pos, 2 all neg
waves=[];
w=[];
d=[];
vpRatio=[];

for k=1:length(allWaveforms)
    theWave=allWaveforms(k).m;
    
    %decide whether to flip (make max pos in every case)
    if theWave(peakPos)<0 & normMode==1
        theWave = theWave*-1;
    end
    if theWave(peakPos)>0 & normMode==2
        theWave = theWave*-1;
    end

    % normalize to peak value=1
    %if normMode==1
    %            theWave = theWave ./ max(theWave);
    %end
    %if normMode==2
    %            theWave = theWave ./ abs(min(theWave));
    %end
    
    % normalize area under the curve = 1
    %totArea = sum(abs(theWave));    
    %theWave = theWave./totArea;
    
    %normalize peak-to-peak amp
    
    r = max(theWave)-min(theWave);
    theWave = theWave./r;
    
    waves(k,:) = theWave;
end

%% get metrics for the waveforms
for k=1:size(waves,1) 
   theWave=waves(k,:);
   
   pos1 = find( theWave == min(theWave) );
   pos2 = find( theWave == max(theWave(pos1(1):end)) );
   
   d(k) = pos2(1)-pos1(1);
   
   metricPos(k,:) = [pos1(1) pos2(1)];
   
   
   %half-valley width
   halfValleyAmp = min(theWave)/2;
   
   beforeIntersect = abs(theWave(1:pos1(1))-halfValleyAmp);
   minPoint = find( beforeIntersect == min(beforeIntersect) );

   afterIntersect = abs(theWave(pos1(1)+1:end)-halfValleyAmp);
   minPointAfter = find( afterIntersect == min(afterIntersect) );
   posAfter = pos1(1)+minPointAfter(1);
   
   w(k) = posAfter-minPoint(1);

   %vpRatio(k) = (theWave(pos2(1)))./(theWave(pos1(1)));
   vpRatio(k) = abs(theWave(pos1(1)))./abs(theWave(pos2(1)));
   
   %debug plot
   debugPlot=0;
   if debugPlot
       figure(10);
       subplot(5,5,k);
       plot(1:256,theWave,'b');
       hold on
       
       plot( minPoint(1), theWave(minPoint(1)), 'or');
       
       
       plot( posAfter, theWave( posAfter ), 'or');
   
       plot( pos1(1), theWave(pos1(1)), 'ok');
       plot( pos2(1), theWave(pos2(1)), 'ok');
       %plot(1:length(beforeIntersect), beforeIntersect, 'r');
       hold off
       title(['w=' num2str(w(k)) ' d=' num2str(d(k)) ' vp=' num2str(vpRatio(k))]);
   end
end

d = d/Fs*1000; %in ms
w = w/Fs*1000;