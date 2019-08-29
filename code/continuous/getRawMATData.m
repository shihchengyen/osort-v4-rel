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

% h = matfile(filename);
h = load(filename); % HM edit
h = h.rw.data; % HM edit

% Check if this file has a time offset variable
if isfield(h,'timestartOffset')
    timestartOffset = h.timestartOffset;
else
    timestartOffset = 0;
end
disp(['OSort Mat format load. Timestart offset used is: ' num2str(timestartOffset)]);

if nargin>1

    try
%         timestamps(:,1) = (1:length(h.analogData)).* (1e6/samplingFreq);
    timestamps(:,1) = (0:(h.analogInfo.NumberSamples-1))' ./ h.analogInfo.SampleRate;
    catch
%         timestamps(:,1) = (1:length(h.analogData)); % for standaloneInit
%         timestamps(:,1) = h.analogTime; % HM edit
    timestamps(:,1) = (0:(h.analogInfo.NumberSamples-1))' ./ h.analogInfo.SampleRate;
    end
    
end

timestamps = timestamps + timestartOffset;
timestamps = timestamps * 1e6; % convert Ripple data (in seconds) to Neuralynx convention (microseconds) % HM edit
timestamps = double(timestamps); % HM edit

if loadScaleFact
   %-- this is for blackrock format
   ElectrodesInfo = h.ElectrodesInfo;
   scaleFact = ElectrodesInfo.MaxDigiValue/ElectrodesInfo.MaxAnalogValue;
   dataSamples = dataSamples./double(scaleFact);
end