%
% locates the closest timestamp to a given value in a series of timestamps
% used to find a datapoint of arbitrary time, where the data is only available for certain times.
%
% mode: 1 use closest (default), 2 use closest which is smaller then the
% requested (begin of block)
%
%urut/april12
function [indMin, tFound] = locateClosestTimestamp( timestamps, tToLocate, mode)
if nargin<3
    mode=1;
end

indMin=[];
tFound=[];

d    = abs(timestamps - tToLocate);  % absolute distance to spike

if mode==1
    indMin = find( min(d)==d );   
else
    % find smallest absolute distance which is also negative
    
    dSigned = (timestamps - tToLocate);  % signed distance to spike
    dSigned_negative = dSigned(find(dSigned<=0));
    
    indMin = find( max(dSigned_negative) == dSigned );   
end

if ~isempty(indMin)
    indMin=indMin(1);
    tFound = timestamps(indMin);
end
