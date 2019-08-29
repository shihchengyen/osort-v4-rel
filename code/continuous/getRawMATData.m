%
%read data from mat file; OSort Mat format.
%
%urut/april07
%kaminskij/ 6 october 2013 reads only part of file that is needed
%urut/april 2015 optimize timestamp setting, add load scale factor for blackrock format
%urut/oct16 added timestartOffset variable to OSort Mat format to allow arbitrary start point in time
%
function	[timestamps,dataSamples] = getRawMATData( filename, fromInd, toInd, samplingFreq, loadScaleFact )
if nargin<5
    loadScaleFact=0;
end

h = load(filename); % HM edit
% h = h.rw.data; % HM edit
h = h.rh.data;

% Check if this file has a time offset variable
if isfield(h,'timestartOffset')
    timestartOffset = h.timestartOffset;
else
    timestartOffset = 0;
end
disp(['OSort Mat format load. Timestart offset used is: ' num2str(timestartOffset)]);

dataSamplesFull(:,1) = h.analogData;
dataSamplesFull = double(dataSamplesFull);
timestampsFull(:,1) = (0:(h.analogInfo.NumberSamples-1))' ./ h.analogInfo.SampleRate;
timestampsFull = double(timestampsFull);

% HM Edit
% For generating a power threshold not biased by large noise, use a 2-step
% process 
% 1. Remove all data larger than 3 s.d. from mean
% 2. Generate a power threshold from there

% % Remove all data points larger than 200uV
% ind_noise = abs(dataSamplesFull) > 200;
% dataSamplesFull = dataSamplesFull(~ind_noise);
% timestampsFull = timestampsFull(~ind_noise);

timestamps = timestamps + timestartOffset;
timestamps = timestamps * 1e6; % convert Ripple data (in seconds) to Neuralynx convention (microseconds) % HM edit
% timestamps = double(timestamps); % HM edit

if loadScaleFact
   %-- this is for blackrock format
   ElectrodesInfo = h.ElectrodesInfo;
   scaleFact = ElectrodesInfo.MaxDigiValue/ElectrodesInfo.MaxAnalogValue;
   dataSamples = dataSamples./double(scaleFact);
end

