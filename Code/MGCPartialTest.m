%% Compute partial multiscale graph correlation and its p-value, given three sample data of size n*p.
%% Still experimental
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
%% @return The partial correlation and its p-value
%%
%% @export

function [corr,pval] = MGCPartialTest(X,Y,Z)
n=size(X,1);
eps=0.000001;

[~,AC]=DCor(X,Z);
[~,BC]=DCor(Y,Z);
[~,CC]=DCor(Z,Z);

[A]=DCorInput(X);
[B]=DCorInput(Y);
[C]=DCorInput(Z);
[A,B]=DCorTransform(A,B,'unbiased',0);
[C,~]=DCorTransform(C,A,'unbiased',0);

PA=A-AC/CC*C;
PB=B-BC/CC*C;
localCor=MGCLocalCor(PA,PB);

% [m,n]=size(localCor);
corr=max(max(localCor));
%[corr,~]=MGCSmoothing(localCor,m,n); 
% if m==1 || n==1
%     corr=localCor(end);
%     optimalScale=m*n;
%     R=1;
% else
%     % find the maximal within the significant region, return optimal scale and a connectec region of local correlations R
%     [corr,optimalScale, R]=MGCSmoothing(localCor,m,n); 
%varPA=MGCLocalCorr(PA,PA)/n/(n-3);
%varPB=MGCLocalCorr(PB,PB)/n/(n-3);
%corr=cov./real(sqrt(varPA*varPB)); % normalized
% corr=corr/size(X,1)^2;% un-normalized

% if varPA<=eps || varPB<=eps
%     corr=0; % if either variance is non-positive, set corr to 0
% end
pval=1-chi2cdf(corr*n+1,1); 
%%% To check validity of chi-square test, verify the following (power should be no more than alpha):
% n=100;p=1;power=0;rep=1000;alpha=0.05;
% for i=1:rep
% X=unifrnd(0,1,n,p);
% Y=unifrnd(0,1,n,p);
% Z=unifrnd(0,1,n,p);
% [corr,pval] = DCorPartialTest(X,Y,Z);
% if pval<alpha;
% power=power+1/rep;
% end
% end
