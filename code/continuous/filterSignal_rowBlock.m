%
%
% filter single trials, aligned in a matrix (one each row)
% also returns the envelope estimate (using hilbert) for each trial
%
function [trialsFiltered,hilbPower] = filterSignal_rowBlock( trials, Hd_LFP2)

trialsFiltered=zeros( size(trials,1), size(trials,2) );
hilbPower=trialsFiltered;
for j=1:size(trials,1)
   trialsFiltered(j,:) = filterSignal( Hd_LFP2, trials(j,:) ); 
   
   hilbPower(j,:) = abs(hilbert(trialsFiltered(j,:)));
end