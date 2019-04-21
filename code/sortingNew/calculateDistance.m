%
%distance measure between different spikes.
%
%
function diffs = calculateDistance(baseClusters, to, weights)

  B = to(ones(size(baseClusters,1), 1), :);% kaminkij 03 07 2017 faster than repmat

diffs = ((baseClusters - B).^2)*weights ; %% kaminkij 03 07 2017 faster than repmat

%%%%diffs = ((baseClusters - repmat(to, size(baseClusters,1),1)).^2)*weights ; %/ (size(baseClusters,2));
