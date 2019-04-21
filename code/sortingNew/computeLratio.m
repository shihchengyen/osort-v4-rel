%
%
%
% computes L-ratio (spike quality metric)
% as defined in Schmitzer-Torbert et al 2005, Neuroscience
%
%
% features used depend on version flag.
% 1: First 3 PC scores, 2 energy normalized version (like in original paper)
%
%
%
function [ L_R ,L,IsolDist, features] = computeLratio( rawWaveforms ,assigned,cluNr, version)
if nargin<4
    version=1;
end

idx=find(assigned==cluNr);

switch(version)
    case 1
        [pc,score,latent,tsquare] = pca(rawWaveforms);
        score=score(:,1:3); %first three PCs
        
        features=score;
    case 2
        %energy of each
        [E,peakAmp,area]=getWaveformEnergy( rawWaveforms);
        E2=repmat(E,size(rawWaveforms,2),1)';
        %normalize each waveform by its energy
        waveforms_normalized=[];
        
        waveforms_normalized=rawWaveforms./E2;

        [pc,score,latent,tsquare] = pca( waveforms_normalized );
    
        features(:,1)=E;
        features(:,2)=peakAmp;
        features(:,3)=area;

        features(:,4)=score(:,1);  % first PC 
        features(:,5)=score(:,2);  %  PC 
        features(:,6)=score(:,3);  %  PC 
        features(:,7)=score(:,4);  %  PC 
        features(:,8)=score(:,5);  %  PC 
        
    otherwise
        error('unknown version')
end

if size(features,2)>size(features,1) | size(features,2)>length(idx)
    %cant compute if more features than spikes are in this cluster
    L=nan;
    L_R=nan;
    IsolDist=nan;
else
    [L,L_R]=L_Ratio( features,idx); %Mclust function 
    IsolDist = IsolationDistance(features, idx);
end



function [E,peakAmp,area] = getWaveformEnergy( waveforms)
E=[];
peakAmp=[];
area=[];
dOfUnits=[];
N=size(waveforms,2);  % nr datapoints per waveform
for j=1:size(waveforms)
    E(j) = sqrt(sum(waveforms(j,:).^2)/N);
    
    [~,Im]=max(abs(waveforms(j,:)));
    
    peakAmp(j) = waveforms(j,Im);  % amplitude of peak, regardless of positive/negative
    area(j) = sum(abs(waveforms(j,:)));
    
        
end