%% Fast and powerful test by subsampling that runs in O(n^2 log(n)+ns*n), based on 
%% C. Shen and J. Vogelstein, “Fast and Powerful Testing for Distance-Based Correlations,”, in preparation
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix, or a m*p data matrix for two sample testing.
%% @param opts - input option structure:
%%        metric is a string that specifies which metric to use, including 'euclidean','hsic', and other variants.%% @param optionCenter is a string that specifies how the distance matrices are centered, including 'mgc'(default), 'unbiased', 'biased' and 'mantel'.
%%        max specifies whether to compute all 1d pairwise statistics, and the number of maximum statistics to use.
%%            0 uses the distance matrix of all dimensions, 1 uses the maximum pairwise stat, pq uses the maximum pq statistics averaged.
%%            Currently it only support 0,1,pq, any other integer in (1,pq) defaults to 1.
%%
%% @return A list contains the following output:
%% @return P-value, test statistic, localCor and optimalScale of MGC (which are 0 for dcor and hsic),
%%
%% @export
%%
function  [stat, pval,localCor,optimalScale]=MGCFastTest(X,Y,opts)
% check input
if nargin<3
    opts = struct('metric','euclidean','max',0); % default parameters
end
if isfield(opts,'metric'); optionMetric = opts.metric; else optionMetric = 'euclidean'; end
if isfield(opts,'max'); optionMax = opts.max; else optionMax = 0; end
[X,Y]=checkTest(X,Y); % check whether it is doing two-sample or independence testing
n=size(Y,1);
p=size(X,2)*size(Y,2);
if optionMax>p
    optionMax=p;
end

% compute the test statistic and optimal scale
opts = struct('metric',optionMetric,'center','mgc','fast',1,'max',optionMax);
[stat,localCor,optimalScale]=MGC(X,Y,opts);

% compute the pvalue via chi-square test
switch optionMax
    case 0
        pval=1-chi2cdf(stat*n+1,1);
    case 1
        pval=1-chi2cdf(stat*n+1,1)^p;
    case p
        pval=1-chi2cdf(stat*n*p+p,p);
    otherwise
        pval=0;
        %         k=sum(sum(corrAll>=corr));
        k=optionMax;
        %         if optionMax<(p/2)
        for j=1:k
            pval=pval+nchoosek(p,j-1)*chi2cdf(corr*n+1,1)^(p-j+1)*(1-chi2cdf(corr*n+1,1))^(j-1);
        end
        %         else
        %             for j=1:(p-k+1)
        %                 pval=pval+nchoosek(p,j-1)*chi2cdf(corr*n+1,1)^(j-1)*(1-chi2cdf(corr*n+1,1))^(p-j+1);
        %             end
        %         end
        pval=1-pval;
end