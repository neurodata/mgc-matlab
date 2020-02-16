%% Transform the distance matrices properly, with column-wise ranking if needed.
%%
%% @param X is an n*n distance matrix;
%% @param Y is an n*n distance matrix;
%% @param optionCenter is a string that specifies how the distance matrices are centered, including 'mgc'(default), 'unbiased', 'biased' and 'mantel'.
%% @param optionRank is a string that specifies whether ranking within column is computed or not.
%%
%% @return A list contains the following output:
%% @return A and B are the centered distance matrices;
%% @return RX and RY are the column rank matrices of X and Y respectively.
%% @return EX and EY are the centering vectors of X and Y depending on the centering scheme.
%%
%% @export
%%
function [A,B,RX,RY,EX,EY]=DCorTransform(X,Y,optionCenter,optionRank)
if nargin<3
    optionCenter='mgc'; % default to mgc center
end
if nargin<4 || strcmpi(optionCenter,'rank')
    optionRank=1; % do ranking or not, 0 to no ranking
end

[A,RX,EX]=DistCentering(X,optionCenter,optionRank);
[B,RY,EY]=DistCentering(Y,optionCenter,optionRank);

%% An auxiliary function that properly transforms the distance matrix X
%%
%% @param X is a symmetric distance matrix;
%% @param optionCenter is a string that specifies how the distance matrices are centered.
%% @param optionRank is a string that specifies whether ranking within column is computed or not.
%%
%% @return A list contains the following output:
%% @return A is the centered distance matrices;
%% @return RX is the column rank matrices of X.
%% @return EX is the centering vector of X depending on the centering scheme.
%%
function [A,RX,EX]=DistCentering(X,optionCenter,optionRank)
[n]=size(X,1);
if optionRank~=0
    RX=DistRanks(X); % the column ranks for X
else
    RX=zeros(n,n);
end

% If rank transformation, take X as the rank matrix then do default mgc transform
if strcmpi(optionCenter,'rank')
    X=RX;
end

switch optionCenter
    case 'mgc' % almost unbiased transform utilized by mgc
        EX=repmat(sum(X,1)/(n-1),n,1);
    case 'unbiased' % unbiased dcor transform
        EX=repmat(sum(X,1)/(n-2),n,1)+repmat(sum(X,2)/(n-2),1,n)-sum(sum(X))/(n-1)/(n-2);
    case 'biased' % biased dcor transform
        EX=repmat(mean(X,1),n,1)+repmat(mean(X,2),1,n)-mean(mean(X));
    case 'mantel' % mantel transform to do simple centering
        EX=sum(sum(X))/n/(n-1);
    otherwise % default to mgc transform
        EX=repmat(sum(X,1)/(n-1),n,1);
        %     case 'dcor' % original dcor that is biased
        %         EX=repmat(mean(X,1),n,1)+repmat(mean(X,2),1,n)-mean(mean(X));
        %     case 'dcor' % single centering of dcor
        %         EX=repmat(mean(X,1),n,1); % column centering
        %     case 'mcor' % single centering of mcor
        %         EX=repmat(sum(X,1)/n,n,1);
        %         EX=EX+X/n;
end
A=X-EX;

% The diagonal entries are excluded for unbiased and mgc centering, but not
% excluded for biased and mantel centering.
if strcmpi(optionCenter,'biased')==0 && strcmpi(optionCenter,'mantel')==0
    for j=1:n
        A(j,j)=0;
    end
end

%% An auxiliary function that sorts the entries within each column by ascending order:
%% For ties, the minimum ranking is used,
%% e.g. if there are repeating distance entries, the order is like 1,2,3,3,4,..,n-1.
%%
%% @param dis is a symmetric distance matrix.
%%
%% @return disRank is the column rank matrices of X.
%%
function [disRank]=DistRanks(dis)

[n,m]=size(dis);
disRank=zeros(n,m);
for i=1:m
    [~,~,a]=unique(dis(:,i));
    disRank(:,i)=a;
end