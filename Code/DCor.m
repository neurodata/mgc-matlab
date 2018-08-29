%% Compute popular global distance or kernel based consistent dependency measure.
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix;
%% @param optionMetric is a string that specifies which global correlation to compute, including 'dcor'(default),'hsic', and other variants.
%% @param optionCenter is a string that specifies how the distance matrices are centered, including 'unbiased'(default) and 'biased' and other variants.
%%
%% @return A list contains the following output:
%% @return corr is the distance correlation in [-1,1], while varX and varY are the distance variances in [0,\infty).
%%
%% @export
%%
function [corr,varX,varY] = DCor(X,Y,optionMetric,optionCenter)
if nargin < 3
    optionMetric='euclidean';
end
if nargin < 4
    optionCenter='unbiased';
end

X=DCorInput(X,optionMetric);
Y=DCorInput(Y,optionMetric);
[A,B]=DCorTransform(X,Y,optionCenter,0);

[corr]=GlobalCov(A,B'); % compute corr statistic
[varX]=GlobalCov(A,A'); % compute varX statistic
[varY]=GlobalCov(B,B'); % compute varY statistic
% varX=diag(varX);
% varY=diag(varY);
corr=corr./real(sqrt(varX*varY'));

if varX<=0 || varY<=0
    corr=0;
end

if strcmpi(optionCenter,'simple')==true || strcmpi(optionCenter,'euclidean2')==true
    corr=abs(corr); % use absolute value for mantel or Pearson correlation via simple or Euclidean2 centering
end

function [cov]=GlobalCov(A,B)
cov=sum(sum(A.*B));