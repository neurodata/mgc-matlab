%% Fast and powerful test by subsampling that runs in O(n^2 log(n)+ns*n), based on 
%% C. Shen and J. Vogelstein, “Fast and Powerful Testing for Distance-Based Correlations,”, in preparation
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix;
%% @param ns specifies the number of subsamples. Generally n/ns should be more than 4, and ns should be large than 10.
%% @param optionMethod is a string that specifies the method, currently supporting 'mgc'(default), 'dcor', 'hsic'.
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
oob='unbiased';
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
        statA(s)=DCor(Xs,Ys,optionMethod,oob);
    end
end

sigma=std(statA)/S;
if strcmpi(optionMethod,'mgc')==true
    mu=max(mean(statA),0);
    [stat,localCor,optimalScale]=MGCSampleStat(X,Y);
else
    mu=0;
    stat=DCor(X,Y,optionMethod,oob);
end
pval=1-normcdf(stat,mu,sigma);
% thres=norminv(1-alpha,0,1);
% ConfidenceInterval=[stat-thres*sigma,stat+thres*sigma];
% if stat<0
%     RequiredSize=inf;
% else
%     RequiredSize=ceil(thres*sigma*n/stat/10)*10;
% end