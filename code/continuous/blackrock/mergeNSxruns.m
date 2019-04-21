%
% merge several osort converted NSx blackrock files into one, establish
% offset between blocks
%
function mergeNSxruns( basepath, prefixStr, channelsToConvert, dirsToMerge, outDir )

nrDirsMerge = length(dirsToMerge);

for k=1:length(channelsToConvert)
    channelNr=channelsToConvert(k);
   
    fname = [ prefixStr num2str(channelNr) '.mat' ];
    
    data_all=[];
    MetaTags_all=[];
    ElectrodesInfo_all=[];
    
    for j=1:nrDirsMerge
        
        %load each and merge
        fNameLoad = [basepath dirsToMerge{j} '/' fname];
        disp(['Loading: ' fNameLoad]);
        load(fNameLoad);
        data_all{j}=data;
        MetaTags_all{j}=MetaTags;
        ElectrodesInfo_all{j}=ElectrodesInfo;
    end
    clear data 
    
    % merge
    blocksize=512;
    
    data=[];
    dataOffsets=[];  % how many datapoints already existed before this block was merged in
    for j=1:nrDirsMerge
        dataToAdd = data_all{j};
        
        %pad the data so it fits in a block of 512
        b = mod(length(dataToAdd),512);        
        dataToAdd = [ dataToAdd zeros(1,512-b) ];      
        
        dataOffsets(j) = length(data); % nr of datapoints that were there before this block started. 
        data = [ data dataToAdd ];        
    end
    
    %now export them
    OSortMATFileversion = 'NSxMerged';
    timestartOffset = 0; % first block starts at 0
    save( [basepath '/' outDir '/' prefixStr num2str(channelNr) '.mat'], 'timestartOffset', 'data', 'dataOffsets', 'dirsToMerge', 'MetaTags', 'ElectrodesInfo','ElectrodesInfo_all', 'MetaTags_all', 'OSortMATFileversion', '-v7.3' );   %v7.3 so file supports partial loading (speedup)
    
end
