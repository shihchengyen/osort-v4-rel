%
%helper function for ground-normalization
%maps an individual channel (wire) to which electrode it belongs; this defines the grouping of channels to macro-electrodes.
%this version groups into 8 wires each; modify if your electrode has more/less channels
%
%returns: elNrOfChannel is the electrode Number (1,2,3,4) for the channel channelNr
%
%urut/july13
function [elNrOfChannel,electrodeAssignment] = mapChannelToElectrode(channelNr)

electrodeAssignment = [ 1:64; repmat(1,1,8) repmat(2,1,8) repmat(3,1,8) repmat(4,1,8) repmat(5,1,8) repmat(6,1,8) repmat(7,1,8) repmat(8,1,8) ]; %which wire is on which electrode

channelInd = find( electrodeAssignment(1,:) == channelNr );

if isempty(channelInd)
	warning(['Channel ' num2str(channelNr) ' not found to assign to electrode']);
end

elNrOfChannel = electrodeAssignment(2,channelInd);


