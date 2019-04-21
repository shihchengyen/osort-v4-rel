
%
% Demonstration of how to merge clusters across different runs (with different detection/alignment parameters)
%

basepath = '~/mnt/dataRawEpilepsy/P33CS_032714/sort_NOv2/';

channelNr=58;
sortDirMin='5.068';
sortDirMax='5.062';
sortDirMerge='5.065';

clustersMin=[765, 775, 797 ];
clustersMax=[759,807, 854];

% list clusters here which are clearly miss-aligned to make spikes available to merge into other clusters. this has to be a subset of what is given in clustersMin/Max
nonPriorityClusters=[807 797];  

priorityFile=2;

clustersMinMaxMerge(basepath, channelNr, sortDirMin, sortDirMax, sortDirMerge, clustersMin, clustersMax, priorityFile,nonPriorityClusters);

