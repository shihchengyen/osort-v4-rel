
% use for realigment of raw spikes after sorting using exact
% method
%
% INPUTS
% filepath  : path to Osort output file. 



function PostSortingRealigne(filepath)


load(filepath);


for i=1:length(useNegative)
    
 
    A=newSpikesNegative(find(useNegative(i)==assignedNegative),:);
    B=mean(A);
   
if B(95)>0
realig  = realigneSpikesWhiten(A, [], 1, 0);  
else
   realig = realigneSpikesWhiten(A, [], 2, 0);   
end

newSpikesNegative(find(useNegative(i)==assignedNegative),:)=realig;
end

save([filepath(1:end-4),'r',filepath(end-3:end)]);