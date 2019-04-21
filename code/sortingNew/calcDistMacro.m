%calculate the appropriate distance for the method chosen
%
%used by sortSpikesOnline.m
%
%urut

% HM Edit: added weight factor to thresholdMethod ~= 1
% To constrain calculation of distance to the region surrounding peak
% (50:200th points) 
function D=calcDistMacro(thresholdMethod, baseSpikes, testSpike, weights,Cinv)

if thresholdMethod==1
    D=calculateDistance(baseSpikes, testSpike, weights);
else
    D=calculateDistanceChi2(baseSpikes, testSpike, Cinv,weights);
end

