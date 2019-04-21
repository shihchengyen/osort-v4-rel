%
% convert a blackrock NS6 (continuous) file to matlab files, one per channel.
% this is a pre-processing step to process NSx files in OSort
%
% fnameNSx: name of the NS6 file
% channelsToConvert: list of channels
% fDirOut: where to store the files
% segNr: which data segment (if several). this refers to NS6.Data{x}. Only one can be exported at a time.
%
% beware - running this requires large amounts of memory
%
% Running this requires 'openNSx' and related code, which is provided by blackrock.
%
%urut/april15
function convertNSx_toMat(fnameNSx, channelsToConvert, fDirOut, segNr, prefixStr)
if nargin<3
    segNr=1;
end
if nargin<5
    prefixStr='BL';
end

for k=1:length(channelsToConvert)
    channelNr = channelsToConvert(k);
    
    NS6=[];
    openNSx('report','read',fnameNSx,['e:' num2str(channelNr)],  'p:double' );
    
    % remove structures, because partial loading does not support them. OSort uses partial loading based on matfile
    if iscell(NS6.Data)
        data = NS6.Data{segNr};   % lowercase for compatibility with OSort original format
    else
        %if only one
        data = NS6.Data;   % lowercase for compatibility with OSort original format
    end

    %data = NS6.Data;
    
    MetaTags = NS6.MetaTags;
    RawData = NS6.MetaTags;
    ElectrodesInfo = NS6.ElectrodesInfo;
    
   
    
    save( [fDirOut '/' prefixStr num2str(channelNr) '.mat'], 'data','MetaTags', 'RawData', 'ElectrodesInfo', '-v7.3' );   %v7.3 so file supports partial loading (speedup)
end