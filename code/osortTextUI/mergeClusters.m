%
%script to manually merge clusters
%
%overwriteParams = [] to request interactively. set to [channel Threshold clsToMerge] for automatic merging
%
%
%urut/feb05
%
function mergeClusters( basePath, sortVersion, figsVersion, finalVersion, overwriteParams )
if nargin<1
    basePath='W:\dataRawEpilepsy\P26CS_121112\';
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

%no parameters below
%========
basePathOut=[basePath '/' sortVersion '/'];
outPath = [basePathOut ];

basePathFigs=[basePath '/' figsVersion '/'];
basePathFigsFinal =[basePathFigs finalVersion '/'];

if exist(basePathFigsFinal)==0
    mkdir(basePathFigsFinal);
end
if exist(outPath)==0
    mkdir(outPath);
end

disp(['mergeClusters.m; currently processing ' basePath '. Using:' basePathOut ' and ' basePathFigs]);

continueLoop=true;
while (continueLoop)
    
    if isempty(overwriteParams)
        %interactive prompt
        channel=input('Channel [Axx]:','s');
        threshold = input('Threshold : ','s');
    else
        channel = num2str(overwriteParams(1));
        threshold = num2str(overwriteParams(2));
        continueLoop=false; % only one run in non-interactive mode
    end
    
    %default prefix is 'A' if none is given
    if ~isempty(str2num(channel(1)))
        channel = ['A' channel];
    end
    
    fname = [basePathOut threshold '/' channel '_sorted_new.mat'];
    
    if exist(fname)==2
        outPath1 = [basePathOut threshold '/'];
        
        disp(['loading:' fname]);
        load(fname);
        
        disp( ['current use (n=' num2str(size(newSpikesNegative,1)) '). Cls=' num2str(useNegative(:)') ]);
        
        if isempty(overwriteParams)
            useNegativeNew = str2num( input('Clusters to merge [a,b] : ','s') );
        else
            useNegativeNew = overwriteParams(3:end); 
        end
        
        disp(['new: ' num2str(useNegativeNew)]);
        
        if length(useNegativeNew)>2
            disp('error, can only merge two clusters at any time');
            continue;
        end
        
        %test whether entered cluster numbers are valid
        if length( intersect( useNegativeNew', useNegative) ) < length(useNegativeNew)
            disp(['Invalid cluster nr: ' num2str(useNegativeNew) ' .Canceled. (Already merged?)']);
            if ~isempty(overwriteParams)
                continueLoop=false;
            end
            
            continue;
        end
        
        if isempty(overwriteParams)
            reply = input('Merge OK? [Y|N] ','s');
        else
            reply='Y';
            
        end
        
        if reply=='Y' || reply=='y'
            
            indsToReplace = find( assignedNegative == useNegativeNew(2) );
            if isempty(indsToReplace)
                warning(['Cl ' num2str(useNegativeNew(2) ) ' empty, cant merge']);
            else
                
                assignedNegative(indsToReplace) = useNegativeNew(1);
                
                %remove the element. dont use setdiff,since it changes the order of the elments
                useNegativeTmp=[];
                for i=1:length(useNegative)
                    if useNegative(i)~=useNegativeNew(2)
                        useNegativeTmp = [ useNegativeTmp ;useNegative(i) ];
                    end
                end
                useNegative = useNegativeTmp;
                
                if exist('useNegativeMerged')
                    useNegativeMerged = [ useNegativeMerged useNegativeNew(2) ];
                else
                    useNegativeMerged = useNegativeNew(2);
                end
                
                fnameOut = [outPath1 channel '_sorted_new.mat'];
                disp(['storing: ' fnameOut]);
                
                %make backup file
                movefile(fnameOut, [fnameOut '.orig']);
                
                save(fnameOut, 'useMUA', 'versionCreated', 'noiseTraces','allSpikesNoiseFree','allSpikesCorrFree','newSpikesPositive', 'newSpikesNegative', 'newTimestampsPositive', 'newTimestampsNegative','assignedPositive','assignedNegative', 'usePositive', 'useNegative', 'useNegativeMerged', 'useNegativeMerged','stdEstimateOrig','stdEstimate','paramsUsed','savedTime');
                disp('file stored, need to manually regenerate figures! ');
            end
            
        end
    else
        warning(['file does not exist ' fname]);
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
