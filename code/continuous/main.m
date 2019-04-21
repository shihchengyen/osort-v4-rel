%run std

for i=9:11
    load('C:\Documents and Settings\Administrator\Desktop\ueli\code\continuous\hd.mat');
    basepath='F:\Data\new\TS_032804\';
    basepathTmp='c:\temp\';
    filename=['a' num2str(i)];
    
    if exist([basepath filename '.mat'])==2
        ['loading ' filename]
        load([basepath filename '.mat']);
        
        [allSpikesPositive,allSpikesNegative, allSpikesTimestampsPositive,allSpikesTimestampsNegative]=processRaw(allSamples, Hd);
        
        
        OKspikesNegative = OKspikesNegative*-1;

        %figure(1)
        %plotWaveformsRaw( OKspikesNegative, OKtimestampsNegative, filename );
        %figure(2)
        %plotWaveformsRaw( OKspikesPositive, OKtimestampsPositive, filename );
        figure(3)
        plotRawWithHist([filename], allSpikesPositive, allSpikesTimestampsPositive)
        
        %extraction finished, store to disk
        save([basepath filename '_spikes_unsorted.mat'], 'OKspikesPositive','OKspikesNegative','OKtimestampsPositive','OKtimestampsNegative');
    end
    clear allSamples;pack;
end


%figure(999)
%plotWaveformsRaw( allSpikesPositive, allSpikesTimestampsPositive, filename );
%figure(1000);
%plotWaveformsRaw( allSpikesNegative, allSpikesTimestampsNegative, filename );
%plotSingleSpikes( OKspikes(1:100,:) );
%plotWaveformsRaw(OKspikes, [],'');

