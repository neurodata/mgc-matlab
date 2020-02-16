%% The main function that tests independent between two data sets by MGC and permutation test.
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix, or a m*p data matrix for two sample testing.
%% @param opts - input option structure:
%%        rep specifies the number of replicates to use for the permutation test, typically >=100.
%%        metric is a string that specifies which metric to use, including 'euclidean','hsic', and other variants.
%%        center is a string that specifies how the distance matrices are centered, including 'mgc', 'unbiased', 'biased' and 'mantel'.
%%        fast specifies whether fast stat computation is used or not. 1 means fast is used, 0 means fast is not used. Not implemented for MGC yet.
%%        max specifies whether to compute all 1d pairwise statistics, and the number of maximum statistics to use.
%%            0 uses the distance matrix of all dimensions, 1 uses the maximum pairwise stat, pq uses the maximum pq statistics averaged, etc.
%%
%% @return A list contains the following output:
%% @return P-value and test statistic of MGC;
%% @return All local p-values by double matrix index, all local correlations by double matrix index, and the estimated optimal scales as matrix single indices.
%%
%% Note that one should avoid report positive discovery via minimizing individual p-values of local correlations,
%% unless corrected for multiple testing problem.
%%
%% @export
%%
function  [statMGC,pMGC,localCor,pLocalCor,optimalScale]=MGCPermutationTest(X,Y,opts)
% check input
if nargin < 3
    opts = struct('rep',100,'metric','euclidean','center','mgc','fast',1,'max',0); % default parameters
end
if isfield(opts,'rep'); rep = opts.rep; else rep = 100; end
if isfield(opts,'metric'); optionMetric = opts.metric; else optionMetric = 'euclidean'; end
if isfield(opts,'center'); optionCenter = opts.center; else optionCenter = 'mgc'; end
if isfield(opts,'fast'); optionFast = opts.fast; else optionFast = 1; end
if isfield(opts,'max'); optionMax = opts.max; else optionMax = 0; end
opts = struct('metric',optionMetric,'center',optionCenter,'fast',optionFast,'max',optionMax);
[X,Y]=checkTest(X,Y); % check whether it is doing two-sample or independence testing
% X=DCorInput(X,option1);
% Y=DCorInput(Y,option1);
np=size(Y,1);

% Compute sample MGC, all local correlations, and the estimated optimal scale
[statMGC,localCor, optimalScale]=MGC(X,Y,opts);
[m,n]=size(localCor);
pLocalCor=zeros(m,n);pMGC=0;

% Compute sample MGC and all local correlations for each permuted data
for r=1:rep
    % Use random permutations on the second data set
    per=randperm(np);
    % Use block permutations
    % per=randpermBlock(np);
    
    YN=Y(per,:);
    [tmp2,tmp]=MGC(X,YN,opts);
    pMGC=pMGC+(tmp2>=statMGC)/rep;
    pLocalCor=pLocalCor+(tmp>=localCor)/rep;
end