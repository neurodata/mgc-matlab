%% Fast and powerful test by subsampling that runs in O(n^2 log(n)+ns*n), based on 
%% C. Shen and J. Vogelstein, “Fast and Powerful Testing for Distance-Based Correlations,”, in preparation
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix;
%% @param ns specifies the number of subsamples. Generally n/ns should be more than 4, and ns should be large than 10.
%% @param optionMethod is a string that specifies the method, currently supporting 'mgc'(default), 'dcor', 'hsic'.
%%% @param optionCenter is a string that specifies how the distance matrices are centered, currently supporting 'unbiased'(default) and 'biased'
%%
%% @return A list contains the following output:
%% @return P-value, test statistic, localCor and optimalScale of MGC (which are 0 for dcor and hsic),
%%
%% @export
%%
%% In comparison, FastTestDCor.m runs in O(ns*n) and the original permutation test runs in O(r n^2 log(n)).
%% This script is slower tha the former, but has a much better power that is similar to the original test.
%%
function  [pval,stat,localCor,optimalScale]=FastTestMGC(X,Y,ns,optionMethod)
% Example 1: ns=20, optionMethod='mgc'; %fast and powerful mgc
% Example 2: ns=20, optionMethod='hsic'; %fast and powerful hsic
% Example 2: ns=20, optionMethod='dcor'; %fast and powerful dcor

if nargin<3
    ns=10; % each subsample has 10 observations by default
end
if nargin<4
    optionMethod='mgc';  % use mgc by default
end
% if nargin<5
optionCenter='unbiased';
% end
n=size(Y,1);
S=floor(n/ns);
% if full data size is no more than 4 times of nn, split to 4 samples --  too few 
% samples will fail the normal approximation and cause the test to be invalid
if n<4*ns 
    ns=floor(n/4);
    S=4;
end
statA=zeros(1,S); % the observed statistics by subsampling

% auxiliary variables and additional processing for MGC subsampling
localCor=0;
optimalScale=0;

% subsampling computation
YP=Y(randperm(n),:);
for s=1:S
    sub=ns*(s-1)+1:ns*s;
    Xs=X(sub,:);
    Ys=YP(sub,:);
    if strcmpi(optionMethod,'mgc')==true
        [statA(s),~]=MGCSampleStat(Xs,Ys);
    else
        statA(s)=DCor(Xs,Ys,optionMethod,optionCenter);
    end
end

% approximate the null distribution by normal distribution
sigma=std(statA)/S;
mu=0;
if strcmpi(optionMethod,'mgc')==true %|| strcmpi(optionCenter,'unbiased')==false
    mu=max(mean(statA),0);
end
% compute the observed statistic
if strcmpi(optionMethod,'mgc')==true
    [stat,localCor,optimalScale]=MGCSampleStat(X,Y);
else
    stat=DCor(X,Y,optionMethod,optionCenter);
end
% compute p value
pval=1-normcdf(stat,mu,sigma);
% if pval<0.05
%     stat
%     mu
%     sigma
% end
% thres=norminv(1-alpha,0,1);
% ConfidenceInterval=[stat-thres*sigma,stat+thres*sigma];
% if stat<0
%     RequiredSize=inf;
% else
%     RequiredSize=ceil(thres*sigma*n/stat/10)*10;
% end