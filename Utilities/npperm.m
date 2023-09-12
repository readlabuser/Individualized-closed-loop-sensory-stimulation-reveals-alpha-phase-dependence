function [p, Ase, Bse] = npperm(A, B, iter, correction, circ, trial_level)

% Last Updated by Tylor J Harlow 5/1/2023.  Analytic assistance from Matthew B Jane
% and Scott C Bressler for benhoch function.
%
% This function performs non-parametric permutation testing of
% time series data. The function is designed to operate at two
% potential levels: individual or group.  
%
% Methodology is based off of: 
%
%   Analyzing Neural Time Series - Cohen 2014
%   Maris & Oostenveld 2007
% 
% Inputs:
%     A & B : Condition arrays of format [subject x time x electrode/frequency].
% 
%     Iter: number of permutation iterations desired. Performed equally at
%     the level of the individual anf the group.
% 
%     Correction: String type.  'benhoch' performs benjamin hochberg
%     correctionn. 'cluster' performs a cluster
%     based correction.
%   
%     Circ: Whether or not to calculate standard error of fischer
%     z-transformation.  Defaults to 0. 
%
%     Trial_level: Whether or not random permutation of trial-level data is
%     to be performed. 


% Setup
N = size(A,1);
N2 = size(B,1);

% Reformat for consistency
dummy = 0;
if ndims(A) < 3
    A = permute(A, [1 3 2]);
    B = permute(B, [1 3 2]);
    dummy = 1;
end

% Convert if data is circular 
if circ == 1
    A = zt(A);
    B = zt(B);
end

% Shuffle Trial Labels at the individual subject level
for j = 1:iter
    if trial_level == 1
        % Shuffle Subjects Sampling From
        idx0 = datasample([1:N]',N,1,"Replace",true);
        idx1 = datasample([1:N2]',N2,1,"Replace",true);
        gA = A(idx0,:,:); gB = B(idx1,:,:);
        grand = cat(1, gA, gB);
        
        % Swap Trial Labels (assume two conditions)
        idx2 = datasample([0:1]',N + N2,1,"Replace",true);
        groupA = grand(find(idx2 == 1),:,:);
        groupB = grand(find(idx2 ~= 1),:,:);

        % Establish null t-distribution
        diff = nanmean(groupA) - nanmean(groupB);
        sd = sqrt((N - 1)*nanstd(groupA) + (N2 - 1)*nanstd(groupB));
        se = sd*sqrt(1/N + 1/N2);
        tstat(j,:,:) = diff./se;
    else
        % Shuffle Subjects Sampling From
        idx0 = datasample([1:N]',N,1,"Replace",true);
        gA = A(idx0,:,:);
        gB = B(idx0,:,:);
        
        % Swap Trial Labels (assume two conditions)
        idx1 = datasample([0:1]',N,1,"Replace",true);
        groupA = cat(1,gA(find(idx1 == 1),:,:),gB(find(idx1 ~= 1),:,:));
        groupB = cat(1,gB(find(idx1 == 1),:,:),gA(find(idx1 ~= 1),:,:));

        % Establish null t-distribution
        diff = groupA - groupB;
        tstat(j,:,:) = nanmean(diff)./(nanstd(diff)./sqrt(N)); 
    end

    % For Cluster Correction
    voxel_pval = 0.05;
    cluster_pval = 0.05;
    tmap = squeeze(tstat(j,:,:));
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % note that here, clusters were obtained by parametrically thresholding
    % the t-maps
    tmap(abs(tmap)<tinv(1-voxel_pval,N-1))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(tmap);
    max_clust_info(j) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps

    % Standard Errors
    Ase(j,:,:) = nanmean(gA);
    Bse(j,:,:) = nanmean(gB);
end


% If less than 3dims, correct previous
if dummy == 1
    A = squeeze(A);
    B = squeeze(B);
end

% Correct Circular Conversion if needed
if circ == 1
    A = inzt(A);
    B = inzt(B);
end

% t-stat is different for unequal number of trials
if trial_level == 1
    % True t-statistic
    truDiff = nanmean(A) - nanmean(B);
    sd = sqrt((N - 1)*nanstd(A) + (N2 - 1)*nanstd(B));
    se = sd*sqrt(1/N + 1/N2);
    truT = truDiff./se;
else
    % True t-statistic
    truDiff = A - B;
    truT = nanmean(truDiff)./(nanstd(truDiff)./sqrt(N));
end



% Signifigance Testing 
t0 = nanmean(tstat);
D = squeeze(truT - t0);
sd = squeeze(nanstd(tstat));
z = (D./sd);  p = 2*normcdf(-abs(z)); 


% Standard Errors
% Condition A
x = squeeze(nanmean(A));
y = squeeze(nanstd(Ase)); 
Aupper = x + y;
Alower = x - y;
% Condition B
x = squeeze(nanmean(B));
y = squeeze(nanstd(Bse));
Bupper = x + y;
Blower = x - y;
if dummy == 1
    Ase = squeeze([Alower'; flipud(Aupper')]);
    Bse = squeeze([Blower'; flipud(Bupper')]);
    p = p';
else
    Ase = squeeze([Alower; flipud(Aupper)]);
    Bse = squeeze([Blower; flipud(Bupper)]);
end


% Cluster-based Correction
if correction == 'cluster'
    zmap = z;
    zmap(abs(zmap) < norminv(1 - voxel_pval)) = 0;
    clustinfo = bwconncomp(z);
    clust_info = cellfun(@numel, clustinfo.PixelIdxList);
    clust_threshold = prctile(max_clust_info, 100 - cluster_pval*100);

    % identify clusters to remove
    whichclusters2remove = find(clust_info<clust_threshold);
    
    % remove clusters
    for i=1:length(whichclusters2remove)
        zmap(clustinfo.PixelIdxList{whichclusters2remove(i)}) = 0;
    end
    p = 2*normcdf(-abs(zmap));  % set p to only signifigant clusters of truT
end



% Benjamin -Hochberg Correction   
if correction == 'benhoch'
    for i = 1:size(p,2)
        p1 = p(:,i);
        [~, CV(i)] = benhoch(p1);
        factor = 0.05/CV(i);
        p(:,i) = p1.*factor;
    end
end



%=================== Benhoch ========================================%
    function [h,CV] = benhoch(p,FDR)
    % [h,CV] = benhoch(p,FDR)
    %   Benjamini-Hochberg Procedure to address false discovery rate (FDR)
    %
    % INPUT VARIABLES
    %     p : vector of individual p-values
    %   FDR : false discovery rate (default = 0.05)
    %
    % OUTPUT VARIABLE
    %   h : Benjamini-Hochberg adjusted hypothesis test
    %       [0=failure to reject null hypothesis, 1=reject null hypothesis]
    %  CV : critical (p) value at which to reject null hypothesis
    %
    % Created: 2018-Sep-07 SCB
    % Revised: 2019-May-17 SCB addressed bug in hypothesis assignment
    
    if nargin<2
        FDR = 0.05; % default FDR value
    end
    
    [pI,I] = sort(p(:)); % sort p-values
    m = numel(pI); % total number of tests
    r = (1:m)'; % inidividual p-value ranks
    CVs = (r/m)*FDR; % Benjamini-Hochberg critical value
    h0 = pI<CVs; % find highest p-value < Benjamini-Hochberg critical value
    h0(1:find(h0==1,1,'last')) = 1; % reset all higher ranks to h = 1
    CV = CVs(find(h0==1,1,'last'));
    if isempty(CV)
        CV = 0.05;
    end
    h(I) = h0; % re-assign hypothesis test results
    h = reshape(h,size(p)); % reshape vector to match original input 'p'
    
    end

    % Fisher z-transform
    function data = zt(data)
        data = (1/2) * log((1 + data) ./ (1 - data)); 
    end

    % Inverse fuscher z-transform
    function data = inzt(data)
        data = (exp(2.*data) - 1)./(exp(2.*data) + 1); 
    end
end
