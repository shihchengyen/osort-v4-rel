%
%produces a figure of the raw data traces and saves it in a file
%
%urut/feb05
%
% using produceRawTraceFig to plot short raw singal with spikes
% kaminskij /feb2014
function  produceRawTraceFigSHORT( handles , figPath, outputFormat,figNR)
if exist(figPath)==0
    mkdir(figPath);
end
if nargin<4
    %default how many plots to make for each cluster
    figNR=1;
end

%get timestamps and header info
handles = initFilter(handles);
[timestampsRaw, nrBlocks,nrSamples,sampleFreq,isContinous,headerInfo, handles.blocksize] = getRawTimestamps( handles.rawFilename, handles.rawFileVersion);
handles.nrSamples=nrSamples;


if length(headerInfo)>0
    %tmp=headerInfo(15);
    %ADbitVolts = str2num(tmp{1}(14:end));
    
    ADbitVolts = str2num(getNumFromCSCHeader(headerInfo, 'ADBitVolts'));
else
    ADbitVolts = 1;
end

%display some statistics about this file
headerInfo

if sampleFreq ~= handles.samplingFreq
    warning(['Sampling rate missmatch between parameters specified and what the file header specifies. This might produce problems, CHECK! Reported: ' num2str(sampleFreq) ' specified ' num2str(handles.samplingFreq) ]);
end

[nrRunsTotal] = determineBlocks( nrSamples );

if min(handles.blockNrRawFig)>nrRunsTotal
    warning('requested block not available - file shorter than block nr');
end

%read raw data
%tillBlocks=handles.blockNrRawFig;

% a=1;
% b=0;
% while b(a)<timestampsRaw(end)
%     a=a+1;
%     b(a)=timestampsRaw(1)+1e6/handles.samplingFreq*handles.blocksize*(a-1);
%  %%%%% b(a)=timestampsRaw(1)+handles.blocksize*(a-1);
%  
% end
%%%%%b=b*handles.samplingFreq
%b(1)=[];
%% find time spike with spike

for i=1:length(handles.useNegative)
    for currentFIG=1:figNR
        A=handles.newSpikesTimestampsNegative( find(handles.assignedClusterNegative==handles.useNegative(i)));
        
        % HM Edit - don't take any spikes that are outside 3 s.d. of
        % cluster
        
        
        ForPLOT(1)=A(3);
        ForPLOT(2)=A(round(length(A)/5));
        ForPLOT(3)=A(round(length(A)/3));
        ForPLOT(4)=A(round(length(A)/2));
        ForPLOT(5)=A(round(length(A)/1.2));
        ForPLOT(6)=A(round(length(A)/3.5));
        ForPLOT(7)=A(round(length(A)/4));
        
        D=find(diff(A)/1e6<0.1);
        if currentFIG~=1
            D(1:(currentFIG-1)*5) =[];
        end
        if length(D)<5
            Endploting=1 ;
        else
            Endploting=0 ;
        end
        ForPLOT(1:length(D))=datasample(A(D),length(D),'Replace',false);
        
        %locate in which block these spikes are located
        for j=1:7
            
            [indMin, tFound] = locateClosestTimestamp( timestampsRaw, ForPLOT(j), 2);
            
            if ~isempty(indMin)
                %map this time to block nrs
                % 1000*512 datapoints are in each block; for each entry in timestampsRaw there are 512 entries; 
                % thus,each block will cover 1000 entries of timestampsRaw
                if handles.rawFileVersion == 5
%                     block(j) = floor(indMin / (512*1000)) + 1;
                    block(j) = floor(indMin / handles.blocksize) + 1;
                else
                    block(j) = floor(indMin / 1000)+1;
                end
            end
            
            %D1=b - ForPLOT(j);
            %D1(D1<0)=NaN;
            %block(j)= find(D1==min(D1));
        end
        
        %%
        figure(888);
        set(gcf,'visible','off');
        close(gcf);   % make sure its closed if it already exists
        figure(888)
        
        set(gcf,'visible','off');
        
        for k=1:length(block)
            
            paramsRaw.howManyBlocks = block(k);
            paramsRaw.startWithBlock = block(k); %only need to read this one block,so start there.
            paramsRaw.includeRange=[];
            %since we arent interested in spikes,following 2 params dont matter for this task
            paramsRaw.prewhiten = 0;
            paramsRaw.alignMethod=1;
            
            paramsRaw = copyStructFields(handles,paramsRaw,{{'paramExtractionThreshold','extractionThreshold'}, ...
                'doGroundNormalization','normalizationChannels', 'pathRaw', 'limit', 'samplingFreq', 'rawFileVersion', 'detectionMethod', 'detectionParams', ...
                'peakAlignMethod','blocksize'});
            
            if isfield(handles,'rawFilePrefix')
                paramsRaw.prefix=handles.rawFilePrefix;
            end
            paramsRaw.rawFilePostfix = handles.rawFilePostfix;
            [allSpikes, allSpikesNoiseFree, allSpikesCorrFree, allSpikesTimestamps, dataSamplesRaw,filteredSignal, rawMean,rawTraceSpikes,runStd2,upperlim, upperlims,stdEstimates, blocksProcessed, noiseTraces, dataSamplesRawUncorrected, blockOffsets, allSpikeInds ] = processRaw(handles.rawFilename, handles.nrSamples, handles.Hd, paramsRaw);
            
            
            
            timeplot=(ForPLOT(k) - blockOffsets)/1e6;
            
            %scale if ADbitVolts value is available
            if ADbitVolts~=1
                filteredSignal=filteredSignal*ADbitVolts*1e6; %convert to uV
                dataSamplesRaw = dataSamplesRaw *ADbitVolts*1e6;
                rawTraceSpikes = rawTraceSpikes * ADbitVolts*1e6;
            end
            
            
            % filter raw signal for illustration purposes (if it is full-band). only enable for debugging
            
            enableBroadbandFilt = 0;
            if enableBroadbandFilt
                n = 4;
                Wn = [2]/(handles.samplingFreq/2);
                [b1,a] = butter(n,Wn,'high');
                HdNew=[];
                HdNew{1}=b1;
                HdNew{2}=a;
                filteredSignal_broadband = filterSignal( HdNew, dataSamplesRaw );
                dataSamplesRaw = filteredSignal_broadband;  %use the filtered version instead (for display purposes only)
            end
            
            
            
            if ~isfield(handles,'plotLimit')
                handles.plotLimit=[1 10]; %in sec
            end
            
            
            plabel=[handles.prefix handles.from ' B' num2str(block(k)) ' T=' num2str(handles.paramExtractionThreshold) ' Fs=' num2str(sampleFreq)];
            try
                try
                    plotSpikeExtractionSHORT(plabel,dataSamplesRaw, dataSamplesRawUncorrected, rawMean, filteredSignal, rawTraceSpikes, runStd2, upperlim, [timeplot-0.01 timeplot+0.09], paramsRaw.samplingFreq, handles, allSpikesTimestamps, blockOffsets, ADbitVolts,k,handles.useNegative(i));
                catch
                    try
                        plotSpikeExtractionSHORT(plabel,dataSamplesRaw, dataSamplesRawUncorrected, rawMean, filteredSignal, rawTraceSpikes, runStd2, upperlim, [timeplot-0.0999 timeplot+0.0001], paramsRaw.samplingFreq, handles, allSpikesTimestamps, blockOffsets, ADbitVolts,k,handles.useNegative(i));
                    catch
                        plotSpikeExtractionSHORT(plabel,dataSamplesRaw, dataSamplesRawUncorrected, rawMean, filteredSignal, rawTraceSpikes, runStd2, upperlim, [timeplot-0.0001 timeplot+0.0999], paramsRaw.samplingFreq, handles, allSpikesTimestamps, blockOffsets, ADbitVolts,k,handles.useNegative(i));
                        
                    end
                end
            end
            
            hold on
            %     plotSpikeExtractionSHORT(plabel,dataSamplesRaw, dataSamplesRawUncorrected, rawMean, filteredSignal, rawTraceSpikes, runStd2, upperlim, handles.plotLimit, paramsRaw.samplingFreq, handles, allSpikesTimestamps, blockOffsets, ADbitVolts);
            
        end
        
        % title(['ADBitsPerVolt ' num2str(ADbitVolts)]);
        
        %title(['Cl ',handles.useNegative(i)])
        scaleFigure;
        %     set(gcf, 'color', 'blue');
        % set(gca, 'color', 'yellow');
        set(gcf, 'InvertHardCopy', 'off');
        %  set(gcf,'Color', [0.900000011920929 0.900000011920929 0.900000011920929]);
        figNameOut=[figPath handles.prefix handles.from '_CL_',num2str(handles.useNegative(i)),'_zRAW',num2str(currentFIG),'.' outputFormat ];
        try
            if outputFormat=='fig'
                disp(['Export FIG: ' figNameOut]);
                saveas(gcf, figNameOut );
            else
                print(gcf, ['-d' outputFormat], figNameOut );
                %              saveas(gcf, figNameOut ,'tiff');
                %              print(gcf, ['-d' 'tiff'], figNameOut );
            end
        catch
            disp('error occured -- see above');
            keyboard;
        end
        
        if isfield(handles,'displayFigures')
            displayFigures=handles.displayFigures;
        else
            %default is close after export to file
            displayFigures=0;
        end
        
        if ~displayFigures
            close(gcf);
        end
        
        if   Endploting==1
            break;
        end
        
    end
end
