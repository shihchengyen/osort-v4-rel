%
%
%part of projection test
%
function projectedResidual = residualProjection( centerDiff, m, spikes )
projectedResidual = zeros(1,size(spikes,1));

%diff between centers, normalize (unit vector)
%centerDiff = m2-m1;
centerDiff = centerDiff/norm(centerDiff);
%norm(centerDiff)

m= repmat(m,size(spikes,1),1);
spikes=spikes-m;
for i=1:size(spikes,1)
 projectedResidual(i)=(spikes(i,:))* centerDiff';
 
 %used to be:
 %     projectedResidual(i) = dot ( (spikes(i,:)-m), centerDiff ) ;

end