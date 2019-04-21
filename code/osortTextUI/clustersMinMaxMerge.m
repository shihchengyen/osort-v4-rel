%Spike sorting tool for OSort; manual merge of a min/max detected sort,
%using findPeak peak alignment method.
%
%The min/max sort should be run using same parameters (detection method,threshold, etc), except alignment (min/max).
%
%priorityFile: 1 neg, 2 pos; this file will be kept and the other merged
%into it. By default, the clusters in the priorityFile take precedence if a
%spike is assigned to two clusters. This can be overwritten by listing
%clusters in the nonPriorityClusters list (see below).
%
%nonPriorityClusters: list of cluster(s) that spikes should be taken away from if there is a conflict.
%this overwrites priorityFile
%
%after running this function, re-run the figures on the resulting file to
%judge the result.
%
%
%urut/may13/CSMC
function clustersMinMaxMerge(basepath, channelNr, sortDirMin, sortDirMax, sortDirMerge, clustersMin, clustersMax, priorityFile, nonPriorityClusters )

% merges a separate min/max merge
% identify clusters to use from min/max that have priority.
% if a spike is assigned to both a min and a max cluster that is listed, it
% will be assigned to the cluster that has priority
%

channelStr = ['A' num2str(channelNr)];
fname1 = [basepath sortDirMin '/' channelStr '_sorted_new.mat'];
fname2 = [basepath sortDirMax '/' channelStr '_sorted_new.mat'];

outDir = [basepath sortDirMerge];
fnameOut = [outDir '/' channelStr '_sorted_new.mat'];
if ~exist(outDir)
    mkdir(outDir);
end

switch(priorityFile)
    case 1
        fnameFirst = fname1;
        fnameSecond = fname2;
        clustersFirst = clustersMin;
        clustersSecond = clustersMax;
    case 2
        fnameFirst = fname2;
        fnameSecond = fname1;
        clustersFirst = clustersMax;
        clustersSecond = clustersMin;
end

if ~exist(fname1) || ~exist(fname2)
    disp(['files dont exist,abort. ' fname1 ' ' fname2]);
end

%first,load the priority file
[spikes1, assigned1, timestamps1, use1, spikesCorrFree1,scalingFactor,noiseTraces] = import_sortResult(fnameFirst);
[spikes2, assigned2, timestamps2, use2, spikesCorrFree2] = import_sortResult(fnameSecond);

[useMUA, versionCreated, stdEstimateOrig,stdEstimate,paramsUsed,savedTime] = import_sortResult_auxiliaryVars(fnameFirst);

% get spikes selected by each
selectedInds1 = selectClusters( clustersFirst, assigned1);
selectedInds2 = selectClusters( clustersSecond, assigned2);

% eliminate overlap
[indsTake1, indsTake2, nrTot, mergeStats] = detectOveralp( timestamps1, selectedInds1, assigned1, timestamps2, selectedInds2, assigned2,nonPriorityClusters );

%display stats of what was merged
if size(mergeStats,1)>0
    mergeIDs = mergeStats(:,1)*1000+mergeStats(:,2);
    for mergeID=unique(mergeIDs)'
        inds= find(mergeIDs==mergeID);
        disp(['From ' num2str(mergeStats(inds(1),1)) ' to ' num2str(mergeStats(inds(1),2)) ' n=' num2str(length(inds)) ]);
    end
end

%verify that overlap is gone
[~, ~, nrTot2] = detectOveralp( timestamps1, indsTake1, assigned1, timestamps2, indsTake2, assigned2, nonPriorityClusters );
if nrTot2 ~=0
    error('overlap not eliminated?');
end

%TODO: re-assign cluster numbers if they are overlapping
if ~isempty(intersect(clustersFirst,clustersSecond))
    error('IDs overlap, need to re-assign');
end

disp(['Tot spikes assigned:' num2str(length(indsTake1)) ' and ' num2str(length(indsTake2)) ' Total ' num2str(length(indsTake1)+length(indsTake2)) ]);

if length(indsTake1)+length(indsTake2) == 0
    error('no spikes assigned - wrong channel? abort to avoid overwrite');
end

% merge
useNegative = [clustersFirst clustersSecond];
newTimestampsNegative = [ timestamps1(indsTake1) timestamps2(indsTake2) ];
newSpikesNegative = [ spikes1(indsTake1,:); spikes2(indsTake2,:) ];

assignedNegative = [ assigned1(indsTake1) assigned2(indsTake2) ];

allSpikesCorrFree = [ spikesCorrFree1(indsTake1,:); spikesCorrFree2(indsTake2,:) ];

usePositive=[];
newTimestampsPositive=[];
assignedPositive=[];
newSpikesPositive=[];
allSpikesNoiseFree=[];

% save
paramsUsedMerge=[];
paramsUsedMerge = addFieldsToStruct(paramsUsedMerge, fnameFirst, fnameSecond, clustersFirst, clustersSecond);

versionCreated='400_minMaxMerge';
paramsUsedOrig=paramsUsed;
paramsUsed=[];
paramsUsed.paramsUsedOrig = paramsUsed;
paramsUsed.sortDirMin = sortDirMin;
paramsUsed.sortDirMax = sortDirMax;
paramsUsed.clustersMin = clustersMin;
paramsUsed.clustersMax = clustersMax;
paramsUsed.priorityFile = priorityFile;
paramsUsed.nonPriorityClusters = nonPriorityClusters;

if exist(fnameOut)
   disp(['File already exists:' fnameOut '. Making backup copy']);
   if exist([fnameOut '.orig'])
       error(['Backup file already exists, abort: ' fnameOut '.orig']);
   else
    movefile(fnameOut, [fnameOut '.orig']);
   end
end

save(fnameOut,'useNegative','newTimestampsNegative','newSpikesNegative','allSpikesCorrFree','paramsUsedMerge','usePositive','newTimestampsPositive', ...
    'newSpikesPositive','allSpikesNoiseFree','scalingFactor','noiseTraces', 'assignedNegative','assignedPositive',...
    'useMUA', 'versionCreated',  'stdEstimateOrig','stdEstimate','paramsUsed','savedTime');

%======= internal functions

function [useMUA, versionCreated, stdEstimateOrig,stdEstimate,paramsUsed,savedTime] = import_sortResult_auxiliaryVars(fname)
load(fname,'useMUA', 'versionCreated', 'stdEstimateOrig','stdEstimate','paramsUsed','savedTime')


function [newSpikesNegative, assignedNegative, newTimestampsNegative, useNegative, allSpikesCorrFree,scalingFactor,noiseTraces] = import_sortResult(fname)
load(fname,'newSpikesNegative', 'assignedNegative', 'newTimestampsNegative', 'useNegative', 'allSpikesCorrFree','scalingFactor','noiseTraces')
if ~exist('scalingFactor')
    scalingFactor=[];
end


function selectedInds1 = selectClusters( clustersToFind, assignedNegative)
selectedInds1=[];
for k=1:length(clustersToFind)
   indsOfCluster = find( assignedNegative == clustersToFind(k) );
   
   if length(indsOfCluster)==0
      warning(['Cluster #' num2str(clustersToFind(k)) ' has 0 spikes. Check cluster number.']); 
   end
   
   selectedInds1 = [selectedInds1 indsOfCluster];
end

%
%timestamps1/2, assigned1/2 are the original (all); selectedInds1/2 are
%used to determine which entries are used.
%
function [indsTake1,indsTake2, nrTot,mergeStats] = detectOveralp( timestamps1, selectedInds1, assigned1, timestamps2, selectedInds2, assigned2, nonPriorityClusters )
mergeStats=[];

tol = 1 * 1000; %in us  ; if timestamps are different by less then this consider them the same.

t1 = timestamps1(selectedInds1);
t2 = timestamps2(selectedInds2);

indsRemoveAll1 = [];
indsRemoveAll2 = [];
%t1 is fixed, eliminate from t2
for k=1:length(t1)
 
    d = t2 - t1(k);
    absD=abs(d);
    
    if min(absD)<=tol
        indsRemove = find( min(absD) == absD );  % indsRemove is an index into t2
        
        % decode where to remove the spike
        cl2 = assigned2(selectedInds2(indsRemove(1)));  %first spike decides
        cl1 = assigned1(selectedInds1(k));
        
        % default is to remove from 2 (cl1 has priority), unless settings indicate that cl1 is non-priority
        useDefaultRemoval=1; %default is to remove from t2
        isCl2_nonPriority = length(find(nonPriorityClusters==cl2));
        isCl1_nonPriority = length(find(nonPriorityClusters==cl1));
        if isCl1_nonPriority && ~isCl2_nonPriority
            useDefaultRemoval=0;
        end
        
        if useDefaultRemoval
            indsRemoveAll2 = [ indsRemoveAll2 indsRemove ];        
            mergeStats(size(mergeStats,1)+1,:) = [cl2 cl1];  %from to
        else
            indsRemoveAll1 = [ indsRemoveAll1 k ];                    
            mergeStats(size(mergeStats,1)+1,:) = [cl1 cl2];  %from to
        end
    end    
end

%inddsRemoveAll1/2 is which entries of selectedInds1/2 to remove (not which
%values)

indsTake1 = selectedInds1(setdiff(1:length(selectedInds1),indsRemoveAll1));
indsTake2 = selectedInds2(setdiff(1:length(selectedInds2),indsRemoveAll2));

nr1 = length(indsRemoveAll1);
nr2 = length(indsRemoveAll2);

nrTot = nr1+nr2;

disp(['nr overlap spikes eliminated: ' num2str(nr1) '/' num2str(nr2)]);

