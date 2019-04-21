%
% overwrite the header of a CSC file. 
% can not exceed 16kb!
%
% originally taken from NRD2CSC.m
%
%urut/feb13
function writeCSCHeader(rFname, headerStr)

if iscell(headerStr)
    % convert header to string
    tHeader = [];
    for i = 1:length(headerStr)
        tHeader = [tHeader, sprintf('%s\r\n',headerStr{i})];    % need \r\n to make header compatible with what is written by neuralynx
    end
    
    headerStr=tHeader;
end


% only overwrite the header part and not the entire file as above (this is tested in UNIX, tbd verify in windows)
fid=fopen(rFname,'r+');  %append mode
fseek( fid, 0, 'bof'); %to beginning
bytesWritten=fprintf(fid, '%s', headerStr)
fclose(fid);