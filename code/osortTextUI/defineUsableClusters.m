%
% script to manually set the valid clusters that should be used for further analysis. sets the variable useNegative in a XX_new_sorted.mat
%
%overwriteParams = [] to request interactively. set to [channel Threshold clsToUse] for automatic 
% call examples: 
% defineUsableClusters( '/Users/Brooks/data/raw/brooks/P28HMH/', 'sortLFP', 'figsLFP', 'final' )
% defineUsableClusters( 'W:\dataRawEpilepsy\P41HMH\', 'sortLFP', 'figsLFP', 'final' )
%
%
%urut/june13: converted to function to avoid accidental re-use of previous variables
%urut/april14: added non-interactive mode
%
function runSuccess=defineUsableClusters( basePath, sortVersion, figsVersion, finalVersion, overwriteParams )
if nargin<1
    basePath='/Users/Brooks/data/raw/';
end
%subdirs to use
if nargin<2
    sortVersion='sortLFP';
end
if nargin<3
    figsVersion='figsLFP';
end
if nargin<4
    finalVersion='final';
end
if nargin<5
    overwriteParams=[];
end
runSuccess=1;

%====
basePathOut=[basePath '/' sortVersion '/'];
outPath = [basePathOut '/' finalVersion '/'];

basePathFigs=[basePath '/' figsVersion '/'];
basePathFigsFinal =[basePathFigs finalVersion '/'];

if exist(basePathFigsFinal)==0
    mkdir(basePathFigsFinal);
end
if exist(outPath)==0
    mkdir(outPath);
end

disp(['defineUsableClusters.m; currently processing ' basePath ' .Using: ' basePathOut ' and ' basePathFigs]);

varList = {'useMUA', 'versionCreated', 'noiseTraces','allSpikesNoiseFree','allSpikesCorrFree','newSpikesPositive', 'newSpikesNegative', 'newTimestampsPositive', 'newTimestampsNegative','assignedPositive','assignedNegative', 'usePositive', 'useNegative', 'useNegativeExcluded', 'stdEstimateOrig','stdEstimate','paramsUsed','savedTime'};
varListToLoad = setdiff(varList, 'useNegativeExcluded');

continueLoop=true;
while (continueLoop)
    
    if isempty(overwriteParams)
        channel = input('Channel [Axx] :','s');
        threshold = input('Threshold : ','s');
    else
        channel = num2str(overwriteParams(1));
        threshold = num2str(overwriteParams(2));       
        continueLoop=false;
    end
    
    %default prefix is 'A' if none is given
    if ~isempty(str2num(channel(1)))
        channel = ['A' channel];
    end
    
    fname = [basePathOut threshold '/' channel '_sorted_new.mat'];
    
    if exist(fname)==2
        
        clear( varList{:} );
        load(fname, varListToLoad{:});  % make sure no old variables linger in the workspace
        
        disp( ['current use (n=' num2str(size(newSpikesNegative,1)) ') == ' num2str(useNegative(:)') ' on Ch:' channel]);
                
        if isempty(overwriteParams)
            useNegativeNew = str2num( input('Clusters to use [a,b,c,d] : ','s') );        
        else
            useNegativeNew = overwriteParams(3:end);
        end
        disp(['New use: ' num2str(useNegativeNew)]);
        
        %test whether entered cluster numbers are valid
        if length( intersect( useNegativeNew, useNegative) ) < length(useNegativeNew)
            
            clsMissing = setdiff(unique(useNegativeNew),useNegative);
            
            disp('error,invalid cluster nr entered. canceled.');
            disp(['Cls not found are: ' num2str(clsMissing)]);
            runSuccess=0;
            %keyboard;
            continue;
        end
        
        %calc stats
        nrTot=0;
        nrTotNew=0;
        for kk=1:length(useNegative)
            nrTot = nrTot + length(find(assignedNegative==useNegative(kk)));
        end
        for kk=1:length(useNegativeNew)
            nrTotNew = nrTotNew + length(find(assignedNegative==useNegativeNew(kk)));
        end
        nrNonNoise = length(find(assignedNegative<99999999));
        
        percDetected = nrTotNew/length(assignedNegative);
        percAssigned = nrTotNew/nrNonNoise;
        disp(['Tot found:' num2str(length(assignedNegative)) ' Tot non-noise:' num2str(nrNonNoise) ' TotAssignedBefore:' num2str(nrTot) ' TotAssignedNew:' num2str(nrTotNew) ' % of detected:' num2str(percDetected) ' % of assigned:' num2str(percAssigned)]);
        
        if isempty(overwriteParams)
            reply = input('New OK? [Y|N] ','s');
        else
            reply='Y';
        end
        
        if reply=='Y' || reply=='y'
            useNegativeOrig=useNegative;
            useNegative=useNegativeNew;
            
            useNegativeExcluded = setdiff(useNegativeOrig,useNegative);
            
            fname=[outPath channel '_sorted_new.mat'];
            disp(['storing: ' fname]);

            % create dummy variables of these dont exist
            if ~exist('useMUA')
                useMUA = 0;
            end
            if ~exist('versionCreated')
                versionCreated = 0;
            end
            if ~exist('stdEstimateOrig')
                stdEstimateOrig = 0;
            end
            if ~exist('stdEstimate')
                stdEstimate = 0;
            end
            if ~exist('paramsUsed')
                paramsUsed = 0;
            end
            if ~exist('savedTime')
                savedTime = 0;
            end
           
            save(fname, varList{:} );
            %'useMUA', 'versionCreated', 'noiseTraces','allSpikesNoiseFree','allSpikesCorrFree','newSpikesPositive', 'newSpikesNegative', 'newTimestampsPositive', 'newTimestampsNegative','assignedPositive','assignedNegative', 'usePositive', 'useNegative', 'useNegativeExcluded', 'stdEstimateOrig','stdEstimate','paramsUsed','savedTime');
            
            %copy the figs, only of the chosen clusters
            for kk=1:length(useNegative)
                
                fnameFrom=[basePathFigs threshold '/' channel '_CL_' num2str(useNegative(kk)) '_THM_*.png'];
                fnameTo=[basePathFigs finalVersion];
                
                disp(['copying from:' fnameFrom ' to ' fnameTo]);
                try
                    copyfile(fnameFrom, fnameTo);
                catch
                end
            end
            disp(['Finished for : ' fname]);
        end
    else
        disp(['file does not exist ' fname]);
        runSuccess=0;
    end
  
    if ~continueLoop  %abort if flag already set to abort
        continue;
    end
    
    continueReply = input('Continue [Y|N] : ','s');
    if continueReply=='Y' || continueReply=='y'
        continueLoop=true;
    else
        continueLoop=false;
    end
end
