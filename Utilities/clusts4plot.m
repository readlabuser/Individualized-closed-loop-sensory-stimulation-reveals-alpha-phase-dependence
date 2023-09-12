function [x, y] = clusts4plot(p, thresh, y)

% Sets up permutation p-values to plot clusters with patch
% Inputs: p = pvalue vector
%         ylim = two element vector of ylimits

% Orient
[r, c] = size(p);
if r > c
    p = p';
end

% Get indices meeting threshold
sig = find(p < thresh);
if isempty(sig)
    fprintf('No significant clusters');
    x = [];
    return
end

% Segment based off of deriviative of indices
segs = [1 diff(sig, 1)];
idx = find(segs > 1);
% If only one cluster
if isempty(idx)
    x = [min(sig) max(sig)];
    x = [x, fliplr(x)];
    y = repelem(ylim, [2 2]);
end

% Starting value
n = 1;

% Loop through cluster count
for i = 1:length(idx) + 1
    
    if i < length(idx) + 1

        % Get cluster i
        tmp = [sig(n) sig(idx(i) - 1)];
        x(i,:) = [tmp, fliplr(tmp)];
        y = repelem(ylim, [2 2]);

        % Set future start value
        n = idx(i);
    else

        % Get cluster end
        tmp = [sig(n) sig(end)];
        x(i,:) = [tmp, fliplr(tmp)];
        y = repelem(ylim, [2 2]);

    end

end

end