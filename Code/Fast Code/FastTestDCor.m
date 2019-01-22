%% Fast Dcor or Hsic test by subsampling that runs in O(ns*n), based on:
%% Q. Zhang, S. Filippi, A. Gretton, and D. Sejdinovic, “Large-scale kernel methods for independence testing,” 
%% Statistics and Computing, vol. 28, no. 1, pp. 113–130, 2018.
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix;
%% @param ns specifies the number of subsamples. Generally n/ns should be more than 4, and ns should be large than 10. 
%% @param optionMethod is a string that specifies the method, supporting 'dcor' and 'hsic'.
%% @param optionCenter is a string that specifies how the distance matrices are centered, including 'unbiased'(default) and 'biased'.
%%
%% @return P-value and test statistic.
%%
%% @export
%%
%% In comparison, FastTestMGC.m runs in O(n^2 log(n)+ns*n) and the original permutation test runs in O(r n^2 log(n)).
%% This script is the fastest one, but has inferior power than the above two at the same sample size.
%%
function  [pval,stat]=FastTestDCor(X,Y,ns,optionMethod,optionCenter)
% Example 1: ns=20, optionMethod='dcor'; %fastest dcor with power loss
% Example 2: ns=20, optionMethod='hsic'; %fastest hsic with power loss

if nargin<3
    ns=10; % each subsample has 10 observations by default
end
if nargin<4
    optionMethod='dcor';  % use mgc by default
end
if nargin<5
    optionCenter='unbiased';  % use unbiased by default
end

n=size(Y,1);
S=floor(n/ns);
% if full data size is no more than 4 times of nn, split to 4 samples --  too few 
% samples will fail the normal approximation and cause the test to be invalid
if n<4*ns 
    ns=floor(n/4);
    S=4;
end
statA=zeros(1,S); % the observed statistics by subsampling

% subsampling computation
YP=Y;
for s=1:S
    sub=ns*(s-1)+1:ns*s;
    Xs=X(sub,:);
    Ys=YP(sub,:);
    statA(s)=DCor(Xs,Ys,optionMethod,optionCenter);
end
mu=0;
sigma=std(statA)/sqrt(S);
stat=mean(statA);
pval=1-normcdf(stat,mu,sigma);