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

% Constrain to window if specified
if nargin > 1
    try
        dataSamples = dataSamplesFull(1,fromInd:toInd);
    catch
        dataSamples = dataSamplesFull(fromInd:toInd,1);
    end
    timestamps = timeStampsFull(fromInd:toInd,1);
else
    dataSamples = dataSamplesFull;
    timestamps = timestampsFull;
end

% if nargin > 1
%     try
%         dataSamples(:,1) = h.analogData(1, fromInd:toInd); % HM Edit
%     catch
%         dataSamples(:,1) = h.analogData(fromInd:toInd,1); % HM Edit
%     end
%     dataSamples = double(dataSamples);
%     
% %     timestamps(:,1) = 1:length(dataSamples); %%%%% WHAT????
% %     timestamps(:,1) = [fromInd:toInd].*(1e6/samplingFreq);
% %     timestamps(:,1) = h.analogTime(fromInd:toInd); % HM Edit
%     timestamps(:,1) = (0:(h.analogInfo.NumberSamples-1))' ./ h.analogInfo.SampleRate;
%     timestamps = timestamps(fromInd:toInd); % constrain to specified block
% else
%     try
% %         timestamps(:,1) = (1:length(h.analogData)).*(1e6/samplingFreq);
%         timestamps(:,1) = (0:(h.analogInfo.NumberSamples-1))' ./ h.analogInfo.SampleRate;
%     catch
% %         timestamps(:,1) = (1:length(h.analogData)); % for standaloneInit
% %         timestamps(:,1) = h.analogTime; % HM Edit
%         timestamps(:,1) = (0:(h.analogInfo.NumberSamples-1))' ./ h.analogInfo.SampleRate;
%     end
% end


timestamps = timestamps + timestartOffset;
timestamps = timestamps * 1e6; % convert Ripple data (in seconds) to Neuralynx convention (microseconds) % HM edit
% timestamps = double(timestamps); % HM edit

if loadScaleFact
   %-- this is for blackrock format
   ElectrodesInfo = h.ElectrodesInfo;
   scaleFact = ElectrodesInfo.MaxDigiValue/ElectrodesInfo.MaxAnalogValue;
   dataSamples = dataSamples./double(scaleFact);
end