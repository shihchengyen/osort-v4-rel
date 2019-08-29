%
%
%mode: 1=with GUI; 2=no GUI (textmode
%
%thresholdMethod: 1=approximation (std), 2=exact (whitened data)
%
function handles = sortMain( hObject, handles, mode, thresholdMethod  )
defineSortingConstants;

handles.correctionFactorThreshold

stdEstimate = handles.stdEstimateOrig + handles.correctionFactorThreshold * handles.stdEstimateOrig;

if thresholdMethod==1
	Cinv=[];
	transformedSpikes=[];
    sortInput = handles.newSpikesNegative;
else
	Cinv = eye(256);

% 	v = 256 -76;% kaminskij when avoiding testing of edges more in sortSpikesOnline exact method
    v = 101; % HM Edit : since we're sorting only on the 50-150th data points = 101 degrees of freedom. 
%     alpha=0.05;	
	alpha=0.05;	% HM Edit
	thres = chi2inv( 1-alpha,v); % Find a value that exceeds 95% of the samples from a chi-square distribution with 76 degrees of freedom. You would observe values greater than 97.3510 only 5% of the time by chance.
	stdEstimate=thres;
    
    sortInput = handles.allSpikesCorrFree;
    transformedSpikes = handles.allSpikesCorrFree;

%     sortInput = handles.newSpikesNegative; % HM Edit
%     transformedSpikes = handles.newSpikesNegative; % HM Edit
   
end

[assigned, nrAssigned, baseSpikes, baseSpikesID] = sortSpikesOnline( sortInput, stdEstimate, size(handles.newSpikesNegative,1), thresholdMethod, Cinv,transformedSpikes, handles.runningAverageLength );

disp(['std estimates ' num2str([stdEstimate handles.stdEstimateOrig]) ] );

%display top-10 clusters max
nrAssignedDisplay=nrAssigned;
if size(nrAssignedDisplay,1)>10
	nrAssignedDisplay=nrAssignedDisplay(end-10:end,:);
end
for j=1:size(nrAssignedDisplay,1)
	disp(['Assigned Cl# ' num2str(nrAssignedDisplay(j,1)) ' n=' num2str(nrAssignedDisplay(j,2))]);
end

minNrSpikes=handles.minNrSpikes;
cluNoise=nrAssigned(find(nrAssigned(:,2)<=minNrSpikes),1);
for i=1:length(cluNoise)
    assigned(find(assigned==cluNoise(i)))=CLUSTERID_NOISE_CLUSTER;
end

cluUse=nrAssigned(find(nrAssigned(:,2)>minNrSpikes),1);
neg='';
for i=1:length(cluUse)
    neg=[neg ',' num2str(cluUse(i))];
end


if mode==1
    set(handles.editUseNeg,'String',neg);
end

handles.useNegative=cluUse;


handles.assignedClusterNegative=assigned;
handles.nrAssignedPreTrim = nrAssigned; % HM Edit
handles.nrAssigned = nrAssigned(ismember(nrAssigned(:,1),cluUse),:); % HM Edit

if mode==1
	guidata(hObject,handles);
end

