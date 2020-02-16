%% Compute all local correlation coefficients in O(n^2 log n)
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix;
%% @param optionMetric is a string that specifies which metric to use, including 'euclidean'(default),'hsic', and other variants.
%% @param optionCenter is a string that specifies how the distance matrices are centered, including 'mgc'(default), 'unbiased', 'biased' and 'mantel'.
%%
%% @return A list contains the following output:
%% @return corr consists of all local correlations within [-1,1] by double matrix index;
%% @return varX contains all local variances for X; varY contains all local covariances for Y.
%%
%% @export
%%
function [corr,varX,varY] = MGCLocalCor(X,Y,optionMetric,optionCenter)
if nargin < 3
    optionMetric='euclidean'; % use euclidean distance by default
end
if nargin < 4
    optionCenter='mgc'; % use mgc centering by default
end

X=DCorInput(X,optionMetric);
Y=DCorInput(Y,optionMetric);

[A,B,RX,RY]=DCorTransform(X,Y,optionCenter);
[cov]=LocalCov(A,B',RX,RY'); % compute all local covariances
[varX]=LocalCov(A,A',RX,RX'); % compute local variances for first data
[varY]=LocalCov(B,B',RY,RY'); % compute local variances for second data
varX=diag(varX);
varY=diag(varY);
% varX=varX(end);
% varY=varY(end);
corr=cov./real(sqrt(varX*varY'));
corr(corr>1)=1; % avoid computational issue that may cause a few local corr to be negligably larger than 1

% set any local correlation to 0 if any corresponding local variance is no larger than 0
for k=1:length(varX)
    if varX(k)<=0
        corr(k,:)=0;
    end
end
for l=1:length(varY)
    if varY(l)<=0
        corr(:,l)=0;
    end
end

%% An auxiliary function that computes all local correlations simultaneously in O(n^2).
%%
%% @param A is the first centered distance matrix;
%% @param B is the second centered distance matrix;
%% @param RX is the column-ranking matrix of A;
%% @param RY is the column-ranking matrix of B.
%%
%% @return covXY is all local covariances computed iteratively.
%%
function [covXY]=LocalCov(A,B,RX,RY)
n=size(A,1);nX=max(max(RX));nY=max(max(RY));
covXY=zeros(nX,nY); %varX=zeros(1,nX); varY=zeros(1,nY);
EX=zeros(1,nX);EY=zeros(1,nY);

% summing up the entrywise product of A and B based on the ranks, which
% yields the local family of covariance and variances
for j=1:n
    for i=1:n
        a=A(i,j);b=B(i,j);k=RX(i,j);l=RY(i,j);
        covXY(k,l)=covXY(k,l)+a*b;
        EX(k)=EX(k)+a;
        EY(l)=EY(l)+b;
    end
end
for k=1:nX-1
    covXY(k+1,1)=covXY(k,1)+covXY(k+1,1);
    EX(k+1)=EX(k)+EX(k+1);
end
for l=1:nY-1
    covXY(1,l+1)=covXY(1,l)+covXY(1,l+1);
    EY(l+1)=EY(l)+EY(l+1);
end
for l=1:nY-1
    for k=1:nX-1
        covXY(k+1,l+1)=covXY(k+1,l)+covXY(k,l+1)+covXY(k+1,l+1)-covXY(k,l);
    end
end

% normalize the covariance by the variances yields the local correlation
covXY=(covXY-EX'*EY/(n-1)/(n));
