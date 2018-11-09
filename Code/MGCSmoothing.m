%% An auxiliary function that finds the optimal local correlation and optimal scale via thresholding and smoothing.
%% Thresholding finds a connected region R of significant local correlations
%% If area of R is too small, return the last local corr; otherwise take the maximum within R.
%%
%% @param localCor is all local correlations;
%% @param m is the number of rows of localCor;
%% @param n is the number of columns of localCor;
%% @param thres is the number of entries for R to be significant;
%%
%% @return A list contains the following output:
%% @return statMGC is the sample MGC statistic within [-1,1];
%% @return optimalScale the estimated optimal scale in matrix single index.
%% @return R is a binary matrix of size m by n indicating the significant region.
%%
function [statMGC, optimalScale,R]=MGCSmoothing(localCor,m,n,sz,thres)
if nargin<4
    sz=min(m,n)-1;
end
if nargin<5
    thres=min(m,n);
end
R=Thresholding(localCor,m,n,sz); % find a connected region of significant local correlations

statMGC=localCor(end); % default sample mgc to local corr at maximal scale
optimalScale=0; % default the optimal scale to 0
if (norm(R,'fro')~=0)
    % tau=0; % number of adjacent scales to smooth with
    if sum(sum(R))>=thres % proceed only when the region area is sufficiently large
        tmp=max(localCor(R==1));
        [k,l]=find((localCor>=tmp)&(R==1)); % find all scales within R that maximize the local correlation
        if tmp >= statMGC
            statMGC=tmp;
            optimalScale=(l-1)*m+k; % take the scale of maximal stat and change to single index
        end
        %         ln=ceil(tau); % number of adjacent rows to check
        %         km=ceil(tau); % number of adjacent columns to check
        %         for i=1:length(k)
        %             ki=k(i);
        %             li=l(i);
        %
        %             % index of adjacent rows and columns
        %             left=max(2,li-ln);
        %             right=min(n,li+ln);
        %             upper=max(2,ki-km);
        %             down=min(m,ki+km);
        %
        %             tmp1=min(localCor(upper:down,li)); % take minimal correlation at given row and along adjacent columns
        %             tmp2=min(localCor(ki,left:right)); % take minimal correlation at given column and along adjacent rows
        %             tmp=max(tmp1,tmp2); % take the min for sample mgc
        %         end
    end
end

%% An auxiliary function that finds a region of significance in the local correlation map by thresholding.
%%
%% @param localCor is all local correlations;
%% @param m is the number of rows of localCor;
%% @param n is the number of columns of localCor;
%%
%% @return R is a binary matrix of size m and n, with 1's indicating the significant region.
%%
function R=Thresholding(localCor,m,n,sz)
thres=max(localCor(end),0);
opt=1; 
if opt==1 % A threshold is estimated based on normal distribution approximation from Szekely2013
    prt=1-0.02/sz; % percentile to consider as significant
    %thres=sqrt(sz*(sz-3)/2-1); % normal approximation, which is equivalent to beta approximation for n larger than 10
    %thres=icdf('normal',prt,0,1)/thres;
    thres=sz*(sz-3)/4-1/2; % beta approximation
    thres=(betainv(prt,thres,thres))*2-1;
    thres=max(thres,localCor(end)); % take the maximal of threshold and local correlation at the maximal scale
end
if opt==2
    thres=localCor;
    thres=thres(thres<0); % all negative correlations
    thres=5*norm(thres,'fro')/sqrt(length(thres)); % 5 times the standard deviation of negative correlations
    thres=max(thres,localCor(end)); % Use the maximal of paratemetric and non-parametric thresholds
end

% Find the largest connected component of significant correlations
R=(localCor>thres);
CC=bwconncomp(R,4);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
if isempty(idx)==false
    CC=CC.PixelIdxList{idx};
    R=zeros(m,n);
    R(CC)=1;
else
    R=0;
end
% Find all correlations that are larger than the threshold
% [~,~,~,R]=FindLargestRectangles((localCor>=thres), [0 0 1],[2,2]);
% % optimalInd=find(optimalInd==1);
% R=double(R);
% if sum(sum(R))==0
%     R=0;
% end