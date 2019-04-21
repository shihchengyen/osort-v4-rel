%
%merges two clusters and saves an updated file.
%this is the function version of the script mergeClusters.m
%
%
function mergeClustersAutomatic(figureoutpath,cluster,channel1,channel2);
display(['merging clusters: ' num2str(cluster) '...']);

global PATH;
outPath = [PATH];

basePathFigs=[figureoutpath];

if exist(outPath)==0
    mkdir(outPath);
end

fname = [PATH channel1];

load(fname);

useNegativeNew = cluster;

%test whether entered cluster numbers are valid
if length( intersect( useNegativeNew, useNegative) ) < length(useNegativeNew)
    error('error,invalid cluster nr entered. canceled.')
end

l = length(useNegativeNew);
for i = 2:l
    indsToReplace = find( assignedNegative == useNegativeNew(i) );
    assignedNegative(indsToReplace) = useNegativeNew(1);

    %remove the element. dont use setdiff,since it changes the order of the elments
    useNegativeTmp=[];
    for j=1:length(useNegative)  % HM Edit. Original code had index 'i', which messed up merging. changed to index 'j'
        if useNegative(j)~=useNegativeNew(i)
            useNegativeTmp = [ useNegativeTmp ;useNegative(j) ];
        end
    end
    useNegative = useNegativeTmp;
end

if exist('useNegativeMerged')
    useNegativeMerged = [ useNegativeMerged useNegativeNew(2:end) ];
else
    useNegativeMerged = useNegativeNew(2:end);
end

% HM Edit - define output filename
extindex = strfind(channel1,'.mat');
channel1short = channel1(1:extindex-1);
outfilename = [channel1short '_merged.mat'];
if exist('useNegativeExcluded') % HM Edit for additional output if merging is done after definition of usable clusters
    save([outPath outfilename], 'useMUA', 'versionCreated', 'noiseTraces','allSpikesNoiseFree','allSpikesCorrFree','newSpikesPositive', 'newSpikesNegative', 'newTimestampsPositive', 'newTimestampsNegative','assignedPositive','assignedNegative', 'usePositive', 'useNegative', 'useNegativeMerged', 'useNegativeMerged','useNegativeExcluded','stdEstimateOrig','stdEstimate','paramsUsed','savedTime','allSpikeInds','nrAssignedPreTrim','nrAssigned'); % HM Edit added output 'allSpikeInds'
else
    save([outPath outfilename], 'useMUA', 'versionCreated', 'noiseTraces','allSpikesNoiseFree','allSpikesCorrFree','newSpikesPositive', 'newSpikesNegative', 'newTimestampsPositive', 'newTimestampsNegative','assignedPositive','assignedNegative', 'usePositive', 'useNegative', 'useNegativeMerged', 'useNegativeMerged','stdEstimateOrig','stdEstimate','paramsUsed','savedTime','allSpikeInds','nrAssignedPreTrim','nrAssigned'); % HM Edit added output 'allSpikeInds'
end

display('finished merging');