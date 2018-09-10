%% MGC and permutation test by fast subsampling. This can also be used for global dependency measure like dcor and hsic to speed up.
%% Note that trivial amount of noise is added to X and Y to break possible ties in data for MGC.
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix;
%% @param ns specifies the number of subsamples. Generally n/ns should be more than 4, and ns should be large than 30.
%% @param optionSubsample specifies the number of replicates to use for the permutation test:
%%        1 uses subsampled statistics for estimating the null only and compute the observed statistic by full data, and runs in O(n^2+ns*n)
%%        2 uses subsampled statistics for estimating the null and also compute the observed statistic by subsampling, which runs in O(ns*n)
%% @param optionMethod is a string that specifies the method, including 'mgc'(default), 'dcor', 'hsic', 'pearson', and other variants.
%% @param alpha specifies the type 1 error level and is used to derive the confidence interval and estimate required sample size to achieve power 1
%%
%% @return A list contains the following output:
%% @return P-value, test statistic, localCor and optimalScale of MGC (which are 0 for non-MGC measures, and size n*n otherwise),
%%         the confidence interval (1*2 vector) for the population correlation with 1-alpha confidence, and the estimate sample size to have power 1 at level alpha,
%%
%% @export
%%
function  [pval,stat,localCor,optimalScale,ConfidenceInterval,RequiredSize]=MGCFastTest(X,Y,ns,optionSubsample,optionMethod,alpha)

% Example 1: ns=50, optionSubsample=2, optionMethod='mgc'; %fastest mgc in linear time, but has less power than original mgc
% Example 2: ns=50, optionSubsample=1, optionMethod='mgc'; %fast mgc saving permutation test and almost no power loss to original mgc
% Example 3: ns=50, optionSubsample=2, optionMethod='dcor'; %fastest dcor with slight power loss
% Example 4: ns=50, optionSubsample=1, optionMethod='dcor'; %fast dcor with almost no power loss
% Example 5: ns=50, optionSubsample=2, optionMethod='hsic'; %fastest hsic with slight power loss
% Example 6: ns=50, optionSubsample=1, optionMethod='hsic'; %fast hsic with almost no power loss

if nargin<3
    ns=100; % each subsample has 100 observations by default
end
if nargin<4
    optionSubsample=1; % use the more powerful version by default
end
if nargin<5
    optionMethod='mgc';  % use mgc by default
end
if nargin<6
    alpha=0.01;
end
n=size(Y,1);
S=floor(n/ns);
% if full data size is no more than 4 times of nn, split to 4 samples --  too few 
% samples will fail the normal approximation and cause the test to be invalid
if n<4*ns 
    ns=floor(n/4);
    S=4;
end
% if rep>0
%     statN=zeros(1,S*rep); % the empirical null statistic if permutation is used.
%     localCorN=zeros(nn,nn,S*rep); % the empirical null statistic if permutation is used.
% end
statA=zeros(1,S); % the observed statistics by subsampling
localCorA=zeros(ns,ns,S); % the local correlations by subsampling
localCor=0;
optimalScale=0;
% add trivial noise to break any ties.
if strcmpi(optionMethod,'mgc')==true
    X=X+(1e-10)*unifrnd(0,1,n,size(X,2));
    Y=Y+(1e-10)*unifrnd(0,1,n,size(Y,2));
end

for s=1:S
    sub=ns*(s-1)+1:ns*s;
    Xs=X(sub,:);
    Ys=Y(sub,:);
    if strcmpi(optionMethod,'mgc')==true
        [statA(s),localCorA(:,:,s)]=MGCSampleStat(Xs,Ys);
    else
        statA(s)=DCor(Xs,Ys,optionMethod);
    end
    %     if rep>0
    %         for r=1:rep
    %             % Use random permutations on the second data set
    %             per=randperm(nn);
    %             YsN=Ys(per,:);
    %             if strcmpi(optionMethod,'mgc')==true
    %                 tmp=MGCSampleStat(Xs,YsN);
    %             else
    %                 tmp=DCor(Xs,YsN,optionMethod);
    %             end
    %             statN((s-1)*rep+r)=tmp;
    %         end
    %     end
end

if strcmpi(optionMethod,'mgc')==true
    sigma=std(statA)/sqrt(S);
    localCor=mean(localCorA,3);
    [stat,optimalScale]=MGCSmoothing(localCor,ns,ns,ceil(ns/sqrt(S)));
    [k,l]=ind2sub([ns,ns],optimalScale);
    k=ceil((k-1)/(ns-1)*(n-1))+1;
    l=ceil((l-1)/(ns-1)*(n-1))+1;
    optimalScale=(l-1)*n+k;
    if optionSubsample<=1 % the observed statistic uses the full data
        [~,localCor,~]=MGCSampleStat(X,Y);
%         [~,optimalScale]=MGCSmoothing(localCor,n,n);
%         [k,l]=ind2sub([n,n],optimalScale)
        stat=localCor(optimalScale);
        sigma=sigma/sqrt(S);
    end
else
    sigma=std(statA)/sqrt(S);
    if optionSubsample<=1
        stat=DCor(X,Y,optionMethod);
        sigma=sigma/sqrt(S);
    else
        stat=mean(statA);
    end
end
pval=1-normcdf(stat,0,sigma);
thres=norminv(1-alpha,0,1);
ConfidenceInterval=[stat-thres*sigma,stat+thres*sigma];
if stat<0
    RequiredSize=inf;
else
    RequiredSize=ceil(thres*sigma*n/stat/10)*10;
end