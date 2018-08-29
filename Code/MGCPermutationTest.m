%% The main function that tests independent between two data sets by MGC and permutation test.
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix;
%% @param rep specifies the number of replicates to use for the permutation test;
%% @param optionMetric is a string that specifies which global correlation to compute, including 'dcor'(default),'hsic', and other variants.
%% @param optionCenter is a string that specifies how the distance matrices are centered, including 'mgc'(default), 'unbiased', 'biased' and 'simple'.
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
function  [pMGC,statMGC,pLocalCor,localCor,optimalScale]=MGCPermutationTest(X,Y,rep,optionMetric,optionCenter)

if nargin<3
    rep=1000; % use 1000 random permutations by default
end
if nargin < 4
    optionMetric='euclidean'; % use euclidean distance by default
end
if nargin < 5
    optionCenter='mgc'; % use mgc centering by default
end
% X=DCorInput(X,option1);
% Y=DCorInput(Y,option1);
np=size(Y,1);

% Compute sample MGC, all local correlations, and the estimated optimal scale
[statMGC,localCor, optimalScale]=MGCSampleStat(X,Y,optionMetric,optionCenter);
[m,n]=size(localCor);
pLocalCor=zeros(m,n);pMGC=0;

% Compute sample MGC and all local correlations for each permuted data
for r=1:rep
    % Use random permutations on the second data set
    per=randperm(np);
    YN=Y(per,:);
    [tmp2,tmp]=MGCSampleStat(X,YN,optionMetric,optionCenter);
    pMGC=pMGC+(tmp2>=statMGC)/rep;
    pLocalCor=pLocalCor+(tmp>=localCor)/rep;
end