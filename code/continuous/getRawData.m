%
%reads raw continous data 
%calls the appropriate routine for the data format used.
%
%the fileFormat variable is defined in defineFileFormat.m (all possible values).
%
%not all formats require a valid samplingFreq argument - it can be set to 0
%in those cases.
%
%urut/april07
function [timestamps,dataSamples] = getRawData( filename, fromInd, toInd, fileFormat, samplingFreq )
if fromInd<1 | toInd <1 %to avoid crash of mex DLL
    error('fromInd/toInd can not be negative');
end

%neuralynx
if fileFormat<=2 || fileFormat==6 || fileFormat==8
    %neuralynx format is in blocks of 512, so divide
    fromInd=ceil(fromInd/512);
    toInd = toInd/512;
    
	[timestamps,dataSamples] = getRawCSCData( filename, fromInd, toInd );
end

%txt file
if fileFormat==3
	[timestamps,dataSamples] = getRawTXTData( filename, fromInd, toInd, samplingFreq );
end

%bin raw
if fileFormat==4
	[timestamps,dataSamples] = getRawBINData( filename, fromInd, toInd, samplingFreq );
end

%MAT file
if fileFormat==5
    %0--> no auto-scaling (raw units as in file)
	[timestamps,dataSamples] = getRawMATData( filename, fromInd, toInd, samplingFreq, 0 );
end

%Blackrock file
if fileFormat==7
    %1--> scale correctly automatically
	[timestamps,dataSamples] = getRawMATData( filename, fromInd, toInd, samplingFreq, 1 );
end
