%
%display the raw trace to compare low and high-gain version for comparison
%of LFP and spike data.
%
%this file assumes 32556Hz sampling rate!!
%
%only use this for manual judgment of LFP/gamma quality.
%
basepath='/data/';

channelNr=6;

patientID='SM_110906'; %p11
%patientID='FC_092206'; %p10


%patientID='SF_032306'; %p9s1
%patientID='SF_032506'; %p9s3

%patientID='LJ_062107'; %p14s1

%patientID='LJ_062307'; %p14s2
%patientID='GV_101207'; %p16
%local
%basepathLFP=['c:/dataLocalCopy/' patientID '/rawLFP/'];
%basepathSpikes=['c:/dataLocalCopy/' patientID '/rawLFP/'];

%server
%rawSubDir='/rawLFP/';
rawSubDir='/raw/';
basepathLFP=['/data/' patientID rawSubDir];
basepathSpikes=['/data/' patientID rawSubDir];

filenameLFP=[basepathLFP 'A' num2str(channelNr) '.Ncs']; %no low cutoff, low-gain signal.
filenameCSC=[basepathSpikes 'A' num2str(channelNr) '.Ncs'];  %high pass,high-gain signal.



cell=['z:/data/sorted/' patientID '/' 'A' num2str(channelNr) '_cells.mat'];
if exist(cell)
    load(cell);
else
    spikes=[0 0 0];
end
    

%basepath='/data2/FC_092206/raw/';
%filename='A10.Ncs';
%events = getRawTTLs([basepath 'events.nev']);

%Fs=32556;
%downsampleFactor=32;
Fs=25000;
downsampleFactor=25;

blockNr=6800;
fromInd=blockNr;
toInd=(blockNr+4000);

[timestamps,nrBlocks,nrSamples,sampleFreq,isContinous,headerInfo] = getRawCSCTimestamps( [ filenameLFP] );
headerInfo

if sampleFreq~=Fs
    error('wrong sampling rate is used');
end

[timestampsLFP,dataSamples] = getRawCSCData( [ filenameLFP], fromInd, toInd );
[timestampsSpikes,dataSamplesCSC] = getRawCSCData( [ filenameCSC], fromInd, toInd );

timestamps=timestampsLFP;  %careful: timestampsLFP and timestampsSpikes sometimes have slight offsets to each other,in which case this is not strictly true.

timestampsRange=[timestamps(1) timestamps(end)];
indsSpikeInRange = find( spikes(:,3)>timestampsRange(1) & spikes(:,3)<timestampsRange(2) );

tmp=headerInfo(15);
ADbitVolts = str2num(tmp{1}(14:end));

%check for out-of-range
maxVal=max(abs(dataSamples));
if maxVal>=32767
    disp('went out of range');
end
maxVal
%convert to volts as measured
dataSamples=dataSamples*ADbitVolts*1e6;


%convert to uV
%dataSamples=dataSamples*1e6;

nrSeconds=2;


figure(90);
subplot(4,1,1)
plot( [(1:Fs*nrSeconds)/Fs]*1000, dataSamples(1:Fs*nrSeconds) )
ylabel('[uV]');xlabel('[ms]');

title([filenameLFP ' raw signal; from LFP extraction']);
xlim([0 nrSeconds*1000]);

n=4; %order of filter


FsRed=Fs/downsampleFactor;
rawDownsampled = resample(dataSamples, 1, downsampleFactor);

%[b,a] = butter(n,10/(Fs/2),'high');
[b1,a1] = butter(n,[3 10]/(FsRed/2)) ;
[b2,a2] = butter(n,[30 100]/(FsRed/2)) ;
[b3,a3] = butter(2,[1 4]/(FsRed/2)) ;
Y1=filtfilt(b1,a1,rawDownsampled);
Y2=filtfilt(b2,a2,rawDownsampled);
Y3=filtfilt(b3,a3,rawDownsampled);

[bspikes,aspikes] = butter(n,[300 3000]/(Fs/2)) ;
Yspikes=filtfilt(bspikes,aspikes,dataSamples);

%filter
%figure(222);  freqz(b,a,1024,FsRed);
subplot(4,1,2)
inds=round(1:FsRed*nrSeconds);
plot( ([1:nrSeconds*Fs]/Fs)*1000, Yspikes(1:nrSeconds*Fs),'g', inds, Y1(inds),'b',  inds, Y2(inds),'r', inds, Y3(inds),'m'  );
legend('spikes','theta','gamma','delta');
title('sw filtered');
xlim([0 nrSeconds*1000]);
ylabel('[uV]');xlabel('[ms]');

subplot(4,1,3)
plot( [(1:Fs*nrSeconds)/Fs]*1000, Yspikes(1:Fs*nrSeconds) )
title([filenameLFP ' spikes, from LFP extraction']);

%--- spikes (as seen on the high-pass filtered/high gain CSC file)
subplot(4,1,4)
[bspikes,aspikes] = butter(n,[300 3000]/(Fs/2)) ;
Yspikes=filtfilt(bspikes,aspikes,dataSamplesCSC);

plot( [(1:Fs*nrSeconds)/Fs]*1000, Yspikes(1:Fs*nrSeconds) )

spikeTimes = spikes(indsSpikeInRange,3);
spikeTimes = spikeTimes-timestampsRange(1);  %normalize
spikeTimes = spikeTimes/1000; %to ms

yval=max(dataSamplesCSC(1:Fs*nrSeconds));

hold on
if ~isempty(spikeTimes)
    plot( spikeTimes, yval, 'x');
end
hold off
xlim([0 nrSeconds*1000]);
title([filenameLFP ' spikes, from SPIKES extraction']);

%% --- powerspectrum

%apply notch filter
%load('~/code/HdNotch_Fs32556_N4.mat');
%load('~/code/HdNotch_Fs25000_N4.mat');
%dataSamples= filtfilthd(HdNotch, dataSamples);

%downsample first
FsDown=1000;
[P,Q] = designDownsampleFactors( Fs, FsDown );
dataSamplesDownSampled = downsampleRawTrace( dataSamples, P, Q ); %FS now 1000Hz

windowSize=FsDown*2; %2s
noverlap=FsDown; %1s
nfft=1024*10;
[pxx,w]=pwelch( dataSamplesDownSampled, windowSize, noverlap,nfft,FsDown,'oneSided');
%[pxx,w]=pmtm( dataSamplesDownSampled, 4,nfft,FsDown);

figure(100);
subplot(1,2,1)
plot(w,log(pxx),'x-')
xlim([1 10])

title(['pwelch ' filenameLFP]);
xlabel('Hz');
ylabel('Power [dB]');

subplot(1,2,2)
plot(w,log(pxx),'x-')
xlim([0 150])

return;

%%

plabel='';
doPlot=1;

indsToUse=100*FsDown:200*FsDown;
figure(300);
[freqs, pestAll, y, A, alpha, coefs] = computeLFPWaveletPOWERSpectrum(dataSamplesDownSampled(indsToUse), plabel, doPlot);
confLim = computeChi2SigOfWaveletpower( y );

[nrFound, allIndsFound, percentages] = computeOscillatoryEpisodes(coefs, confLim, freqs);

[dataFilteredLowPass1, filterSettings2] = filterLFPofTrial( {dataSamplesDownSampled}, 3, FsDown );
[dataFilteredLowPass2, filterSettings2] = filterLFPofTrial( {dataSamplesDownSampled}, 2, FsDown );

figure(202);
subplot(2,2,1)
plot(freqs,nrFound);
subplot(2,2,2)
plot(freqs,percentages);

figure(200);
subplot(4,1,1)
plot( dataSamplesDownSampled(indsToUse) );
subplot(4,1,2);
plot( dataFilteredLowPass1{1}(indsToUse) );
subplot(4,1,3);
plot( dataFilteredLowPass2{1}(indsToUse) );
subplot(4,1,4);
plot( abs(coefs(9,indsToUse)));

