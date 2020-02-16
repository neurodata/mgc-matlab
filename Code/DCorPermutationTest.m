%% Standard permutation test for independence testing using distance correlation, runs in O(rep * n^2).
%% Generally slow and not recommended unless n<100. Use DCorFastTest.m for large sample size.
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix, or a m*p data matrix for two sample testing.
%% @param opts - input option structure:
%%        rep specifies the number of replicates to use for the permutation test, typically >=100.
%%        metric is a string that specifies which metric to use, including 'euclidean','hsic', and other variants.
%%        center is a string that specifies how the distance matrices are centered, including 'unbiased', 'biased', 'mantel' and 'mgc'.
%%        fast specifies whether fast stat computation is used or not. 1 means fast is used, 0 means fast is not used.
%%        max specifies whether to compute all 1d pairwise statistics, and the number of maximum statistics to use.
%%            0 uses the distance matrix of all dimensions, 1 uses the maximum pairwise stat, p uses the maximum p statistics averaged, etc.
%%
%% @return A list contains the following output:
%% @return P-value and test statistic.
%%
%% @export
%%
function  [corr,pval]=DCorPermutationTest(X,Y,opts)
% check input
if nargin < 3
    opts = struct('rep',100,'metric','euclidean','center','unbiased','fast',1,'max',0); % default parameters
end
if isfield(opts,'rep'); rep = opts.rep; else rep = 100; end
if isfield(opts,'metric'); optionMetric = opts.metric; else optionMetric = 'euclidean'; end
if isfield(opts,'center'); optionCenter = opts.center; else optionCenter = 'unbiased'; end
if isfield(opts,'fast'); optionFast = opts.fast; else optionFast = 1; end
if isfield(opts,'max'); optionMax = opts.max; else optionMax = 0; end
opts = struct('metric',optionMetric,'center',optionCenter,'fast',optionFast,'max',optionMax);
[X,Y]=checkTest(X,Y); % check whether it is doing two-sample or independence testing
n=size(X,1);

% Calculate the observed test statistics for the given data sets
pval=0;
if strcmpi(optionMetric,'hhg')
    ch=1;
else
    ch=2;
end

switch ch
    case 1
        corr=HHG(X,Y);
    case 2
        corr=DCor(X,Y,opts);
end

% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % Use random permutations;
    per=randperm(n);
    % Use block permutations
    % per=randpermBlock(n);
    
    DN=Y(per,:);
    switch ch
        case 1
            tmp=HHG(X,DN);
        case 2
            tmp=DCor(X,DN,opts);
    end
    pval=pval+(tmp>=corr)/rep;
end