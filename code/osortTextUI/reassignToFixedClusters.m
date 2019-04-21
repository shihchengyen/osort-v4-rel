%
%this file is used to re-assign all detected waveforms on a channel to a set of pre-defined cluster centers. no new clusters are created.
%
%input is the reassign structure, which describes the channel and the clusters that serve as new centers.
%
%urut/sept06

function reassignToFixedClusters

reassign(1).channel=47;
reassign(1).thres=5;
reassign(1).thresOut=5.01;

reassign(1).validClusters =[1166 1277 1307 1309 1022 1271 1180];


basedir='~/mnt/dataRawEpilepsy/P29CS_103013/';
basedirIn='sort_NOr';
basedirOut='sort_NOr';


for i=1:length(reassign)
    channelNr = reassign(i).channel;
    thres = reassign(i).thres;
    thresOut=reassign(i).thresOut;
    clus = reassign(i).validClusters;
    
    %open files
    fname=[basedir basedirIn '/' num2str(thres) '/A' num2str(channelNr) '_sorted_new.mat'];
    
    fnameOut=[basedir basedirIn '/' num2str(thresOut) '/A' num2str(channelNr) '_sorted_new.mat'];

    load(fname);
    
    
    %find mean waveforms
    meanClusters=[];
    totN=0;
    for j=1:length(clus)
        indsClus = find( clus(j) == assignedNegative );
        totN=totN+length(indsClus);
        meanClusters(j,:) = mean( newSpikesNegative(indsClus,:) );
    end
    
    figure(1)
    plot( meanClusters');
    
    %now re-assign to those,sorting without creating new clusters
    assigned = assignSpikesToClusters(meanClusters, clus, newSpikesNegative, stdEstimateOrig);
    
    totN2=0;
    for j=1:length(clus)
        totN2 = totN2 + length(find(assigned==clus(j)));
    end

    figure(5);
    plotSortingResultRaw(newSpikesNegative, [], assigned, [], [clus ], [], '', {'r','g','b','y','k','m','c'} )
    
    n=size(newSpikesNegative,1);
    [totN totN2 n totN/n totN2/n]
    
    useNegative=clus;
    assignedNegative=assigned;
    
    save(fnameOut);
end
