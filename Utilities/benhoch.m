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
r = (1:m)'; % individual p-value ranks
CVs = (r/m)*FDR; % Benjamini-Hochberg critical value
h0 = pI<CVs; % find highest p-value < Benjamini-Hochberg critical value
h0(1:find(h0==1,1,'last')) = 1; % reset all higher ranks to h = 1
CV = CVs(find(h0==1,1,'last'));
if isempty(CV)
    CV = min(CVs);
end
h(I) = h0; % re-assign hypothesis test results
h = reshape(h,size(p)); % reshape vector to match original input 'p'
end
