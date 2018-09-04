%% The main function that computes the MGC measure between two datasets:
%% It first computes all local correlations,
%% then use the maximal statistic among all local correlations based on thresholding.
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix;
%% @param optionMetric is a string that specifies which global correlation to compute, including 'dcor'(default),'hsic', and other variants.
%% @param optionCenter is a string that specifies how the distance matrices are centered, including 'mgc'(default), 'unbiased', 'biased' and 'simple'.
%%
%% @return A list contains the following output:
%% @return statMGC is the sample MGC statistic within [-1,1];
%% @return localCor consists of all local correlations by double matrix index;
%% @return optimalScale the estimated optimal scale in matrix single index, 
%%         which can be changed to double indices by [k,l]=ind2sub(size(localCor),optimalScale);
%%
%% @export
%%
function [statMGC, localCor, optimalScale,R]=MGCSampleStat(X,Y,optionMetric,optionCenter)
if nargin < 3
    optionMetric='euclidean'; % use euclidean distance by default
end
if nargin < 4
    optionCenter='mgc'; % use mgc centering by default
end
localCor=MGCLocalCor(X,Y,optionMetric,optionCenter); % compute all localCor
[m,n]=size(localCor);
if m==1 || n==1
    statMGC=localCor(end);
    optimalScale=m*n;
    R=1;
else
    % find the maximal within the significant region, return optimal scale and a connectec region of local correlations R
    [statMGC,optimalScale, R]=MGCSmoothing(localCor,m,n); 
end

