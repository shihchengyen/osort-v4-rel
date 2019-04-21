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

h = matfile(filename);

% Check if this file has a time offset variable
if ~isempty(whos(h,'timestartOffset'))
    timestartOffset = h.timestartOffset;
else
    timestartOffset = 0;
end
disp(['OSort Mat format load. Timestart offset used is: ' num2str(timestartOffset)]);

if nargin>1
    try
    dataSamples(:,1) = h.data( 1,fromInd:toInd );
    catch
        dataSamples(:,1) = h.data( fromInd:toInd ,1);  
    end
    
    %timestamps(:,1) = 1:length(dataSamples);
    timestamps(:,1) =  [fromInd:toInd].*(1e6/samplingFreq);
    
else
    
    try
        timestamps(:,1) = (1:length(h.data)).* (1e6/samplingFreq);
    catch
        timestamps(:,1) = (1:length(h.data)); % for standaloneInit
    end
    
end

timestamps = timestamps + timestartOffset;

if loadScaleFact
   %-- this is for blackrock format
   ElectrodesInfo = h.ElectrodesInfo;
   scaleFact = ElectrodesInfo.MaxDigiValue/ElectrodesInfo.MaxAnalogValue;
   dataSamples = dataSamples./double(scaleFact);
end