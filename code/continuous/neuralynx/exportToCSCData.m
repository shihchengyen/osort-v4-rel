%
%export simulated data to Ncs (neuralynx) files.
%
%automatically reshapes matrices appropriately; 
%
%fname filename
%Fs is sampling rate
%data is a simple-by-sample (continous) of the data
%offset in us (timestamp of the first datapoint)
%
%Data needs to be in the appropriate range (uV); data <1 is rounded to 0.
%
%returns the data as written (normally not of interest, for debugging).
%
%This function only works on windows, due to Windows-specific neuralynx mex
%files.
%
%headerLines: cell array of header lines
%
%urut/nov08
%urut/aug11 updated to mex files v5, correctly write sampling rate as well.
function [timestamps,Samples] = exportToCSCData(fname, Fs, data, offset, channelNr, origInfo, headerLines)
if nargin<7
    headerLines={};
end
    
N=length(data);

blocksize=512;

%shorten data if not multiple of 512
data = data(1:end-rem( N, 512) ); %needs to be a multiple of 512

traceNrSamples=length(data);

dataReshapped = reshape(data, blocksize, traceNrSamples/blocksize);

stepsize=1e6/Fs; %in us
timestamps=[1:stepsize*blocksize:traceNrSamples*stepsize]; % in us
timestamps=timestamps+offset; %some arbitrary offset
timestamps=round(timestamps);

AppendFile=0;
ExtractMode=1;
ModeArray =1;
NumRecs=   length(timestamps); %nr timestamps

%all
FieldSelection(1) = 1; %timestamps
FieldSelection(2) = 0; %channel nr
FieldSelection(3) = 0; %sample freq
FieldSelection(4) = 0; % valid samples
FieldSelection(5) = 1; %samples
FieldSelection(6) = 1;%header

%Timestamps = timestamps;
ChannelNumbers=repmat(channelNr,1,NumRecs);
SampleFrequencies=repmat(Fs,1,NumRecs);
NumberOfValidSamples=repmat(blocksize,1,NumRecs);

Samples = dataReshapped;

clock1=clock;
%processedOn = getProcessedOnString;

%somehow only first line of header is exported, unclear why
header{1} = ['######## sim Data from ' origInfo ' Fs=' num2str(Fs) ' urut-Aug11'];
header{2} = ['## XX'];
header{3} = [''];
%header{4} = [''];
%header(1)='line1';
%header(2)='line2';

%uses the new v5 mex files from neuralynx
if computer=='PCWIN64' | computer=='PCWIN32'
    % Windows version
    Mat2NlxCSC(fname, 0, 1, 1, [1 1 1 1 1 1], timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, header);
else
    %Unix version
    %NumRecs=length(timestamps);
    %Mat2NlxCSC(fname, 0, 1, 1, NumRecs, [1 1 1 1 1 1], timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, header);

    disp(['Using Unix Version of Mat2NlxCSC']);
    putRawCSC( fname, timestamps, data,ChannelNumbers, SampleFrequencies, NumberOfValidSamples, headerLines );

    %Function( Filename, AppendFile, ExtractMode, ModeArray, NumRecs, FieldSelection, Timestamps, ChannelNumbers, SampleFrequency, NumberValidSamples, Samples, Header );
           
end

