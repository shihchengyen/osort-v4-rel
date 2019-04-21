%find all projection test values (distances and fit)
%
%sessions: string identifier of sessions to consider
%considerAreas: id of areas that should be considered. neurons from other areas are excluded.
%
%
%urut/feb05
%urut/march07. revised for recall paper.
%urut/MPI/oct11 added waveform analysis
%urut/CSMC/jan13 added new ISI variability metrics CV,CV2, etc
%
function [statsProj, nrValidChannels, ISIstats,basedirsUsed, allWaveforms] = evalProjectionTestAllCells( sessions, sessionsTaskname, considerAreas, sortSubdir )
nrValidChannels=0;

path='~/mnt/dataPrevious/';
path2='~/mnt/dataPrevious/dataBubbles/';

if nargin<4
    sortSubdir='sort';
end

%spike sorting paper
%sessions={'HM_102104','HM_102204'};
%sessionID=[301 302];  % patient # * 100 + session #

%recall paper

%sessionID, channel, clusterNr1, clusterNr2, distance, fit1, fit2, mode(1: >1 neuron, 2: only 1 neuron, d + fit2 invalid)
statsProj=[];
ISIstats=[];
allWaveforms=[];
for z=1:length(sessions)
    basedir=[path sessions{z} '/' sortSubdir sessionsTaskname{z} '/final/'];
    %basedirRaw=[path sessions{z} '/raw/'];

    basedirsUsed{z} = basedir;

    fNameBrainArea = [path '/events/' sessions{z} '/' sessionsTaskname{z} '/brainArea.mat'];
    if exist(fNameBrainArea)
        load(fNameBrainArea);
    else
        %alternative location for brainarea file
        fNameBrainArea2 = [path sessions{z} '/' sessionsTaskname{z} '/brainArea.mat'];
        if exist(fNameBrainArea2)
            load(fNameBrainArea2);
        else
           %error(['brainArea.mat not found in both: '  fNameBrainArea ' and ' fNameBrainArea2]);
           keyboard;
        end
        
    end
    
    for channel=1:64
        
        fname=[basedir 'A' num2str(channel) '_sorted_new' '.mat'];
        
        %-- test whether file exists (has any neurons)
        if exist(fname)~=2
            disp([fname ' does not exist,skip']);
            continue;
        end

        %-- test whether neuron is in an area that we want to consider
        indsArea = find( brainArea(:,1) == channel);
        
        areaOfChannel = brainArea(indsArea(1),4);
        
        if length ( find( considerAreas == areaOfChannel) ) == 0
            disp(['ignore channel, wrong area. ' sessions{z} ' C:' num2str(channel) ' A:' num2str(areaOfChannel)]);
            continue;
        end

        
        %-- try to load the file
        nrValidChannels=nrValidChannels+1;

        %disp(['would load: ' fname]);        
        %continue;
        
        
        load(fname);   % load newSpikesNegative, assignedNegative  (the spikes)
        disp(['processing : ' fname]);

        %-- ISI
        [statsTmp, nrSorted, allWaveformsOfFile] = getSNRISISortingStatsOfCluster(channel, z, useNegative, assignedNegative, newSpikesNegative, newTimestampsNegative, stdEstimateOrig);
        ISIstats=[ISIstats; statsTmp];
        
        
        allWaveforms = [ allWaveforms allWaveformsOfFile ];
        
        %-- projection test
        clusters = useNegative;
        %find all pairs
        pairs=[];
        c=0;
        %significance test between all clusters
        for i=1:length(clusters)
            for j=i+1:length(clusters)
                c=c+1;
                pairs(c,1:2)=[clusters(i) clusters(j)];
            end
        end

        %if only 1 cluster,find fit in any case
        if length(clusters)==1
            pairs(1,1:2)=[clusters(1) 0];
        end

        for k=1:size(pairs,1)
            clNr1 = pairs(k,1);
            clNr2 = pairs(k,2);
            disp(['Pair: ' num2str([channel clNr1 clNr2]) ]);

            [d,residuals1,residuals2,Rsquare1, Rsquare2] = figureClusterOverlap(allSpikesCorrFree, newSpikesNegative, assignedNegative, clNr1, clNr2, '', 3, '');

            %if only 1 cluster,use only goodness of fit and not distance
            mode=1;
            if clNr2==0 %      size(pairs,1)==1
                mode=2;
                d=0;
                Rsquare2=0;
            end

            entry = [ z channel clNr1 clNr2 d Rsquare1 Rsquare2 mode];
            entryNr = size(statsProj,1)+1;
            statsProj(entryNr,:) = entry;
        end
    end
end

