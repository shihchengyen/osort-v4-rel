function [normalizationChannels,paramsIn,filesToProcess] = StandaloneGUI_prepare(noiseChannels,doGroundNormalization,paramsIn,filesToProcess,filesAlignMax, filesAlignMin, normalizeOnly, groundChannels )

%kaminskij - if runningAverageLength parameter is not provided 
if ~isfield(paramsIn, 'runningAverageLength')
paramsIn.runningAverageLength=100;
end
%kaminskij - if blocksize parameter is not provided 
if ~isfield(paramsIn , 'blocksize')
 paramsIn.blocksize=512000;   
end

normalizationChannels=[]; %maps channel to electrode. All that are listed are used for normalization.
if exist('noiseChannels')
    filesToProcess=setdiff(filesToProcess, noiseChannels);
end

if doGroundNormalization

    %map channels to El number
    for j=1:length(normalizeOnly)
	elNrOfChannel = mapChannelToElectrode(normalizeOnly(j));
        normalizationChannels(1:2,j) = [normalizeOnly(j); elNrOfChannel];
    end
    
    %for j=1:length(filesToProcess)
    %    if length(find(normalizeOnly==filesToProcess(j)))==1
    %        normalizationChannels(1:2,j) = [filesToProcess(j); electrodeAssignment(2, j)];
    %    end
    %end

    excludeChannels=[]; %setdiff(filesToProcess, filesToProcess);

    if exist('noiseChannels')
        excludeChannels=[excludeChannels noiseChannels];
    end

    %exclude noise channels from the grand average.
    excludeInds=[];
    for i=1:length(excludeChannels)
        excludeInds = [ excludeInds find( normalizationChannels(1,:) == excludeChannels(i) ) ];
    end
    includeInds = setdiff( 1:size(normalizationChannels,2), excludeInds);
    normalizationChannels = normalizationChannels(:, includeInds);

    paramsIn.groundChannels=groundChannels;
end

%align method can be changed for each channe;; alignMethod is only used if peakAlignMethod=1 (findPeak)
paramsIn.alignMethod = repmat(paramsIn.defaultAlignMethod, 1, length(filesToProcess));
for i=1:length(filesToProcess)
	if length(find( filesAlignMax == filesToProcess(i) ))==1
		paramsIn.alignMethod(i) = 1;
	end
	if length(find( filesAlignMin == filesToProcess(i) ))==1
		paramsIn.alignMethod(i) = 2;
	end
end
