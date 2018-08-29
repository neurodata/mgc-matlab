%% The function that tests independent between two data sets by distance-based correlation and permutation test.
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix;
%% @param rep specifies the number of replicates to use for the permutation test;
%% @param optionMetric is a string that specifies which global correlation to compute, including 'dcor'(default),'hsic', and other variants.
%% @param optionCenter is a string that specifies how the distance matrices are centered, including 'unbiased'(default) and 'biased' and other variants.
%%
%% @return A list contains the following output:
%% @return P-value and test statistic.
%%
%% @export
%%
function  [pval, test]=DCorPermutationTest(X,Y,rep,optionMetric,optionCenter)
if nargin<3
    rep=1000;
end
if nargin<4
    optionMetric='euclidean';
end
if nargin<5
    optionCenter='unbiased';
end
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
        test=HHG(X,Y);
    case 2
        test=DCor(X,Y,optionMetric,optionCenter);
end

% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % Use random permutations;
    per=randperm(n);
    DN=Y(per,:);
    switch ch
        case 1
            tmp=HHG(X,DN);
        case 2
            tmp=DCor(X,DN,optionMetric,optionCenter);
    end
    pval=pval+(tmp>=test)/rep;
end