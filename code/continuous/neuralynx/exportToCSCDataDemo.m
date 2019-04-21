%demonstrates how to use exportToCSCData.m to export simulated data.
%

traceLength=5;
Fs=32556;
data=round( cumsum(randn(1, Fs*traceLength)) );  

%export data to neuralynx Ncs files
fname='/tmp/exp.Ncs';

offset=500000;
[timestampsW,SamplesW] =exportToCSCData(fname, Fs, data, offset, 5, 'origfname.mat');

%% read back for testing
[timestamps,nrBlocks,nrSamples,sampleFreq,isContinous,headerInfo] = getRawCSCTimestamps( fname );
headerInfo


[timestamps2,dataSamples] = getRawCSCData( fname, 1, 999 );

figure(1);
subplot(1,2,1);
plot( data );
title('simulated');
subplot(1,2,2);
plot(  dataSamples );
title('read back from file');



%%
fname='/tmp/testCSC/t1-CSC1.ncs';
fname='/media/DataH/humanRecordings/p27cs_2013-01-19_11-17-49_ueli/CSC1.ncs';

[timestamps,nrBlocks,nrSamples,sampleFreq,isContinous,headerInfo] = getRawCSCTimestamps( fname );
headerInfo

[timestamps2,dataSamples] = getRawCSCData( fname, 1, 999 );

figure(10);
plot( dataSamples.*  0.006104);