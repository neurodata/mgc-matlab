%% The main function that computes the MGC measure between two datasets:
%% It first computes all local correlations,
%% then use the maximal statistic among all local correlations based on thresholding.
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix, or a m*p data matrix for two sample testing.
%% @param opts - input option structure:
%%        metric is a string that specifies which metric to use, including 'euclidean','hsic', and other variants.
%%        center is a string that specifies how the distance matrices are centered, including 'mgc', 'unbiased', 'biased' and 'mantel'.
%%        fast specifies whether fast stat computation is used or not. 1 means fast is used, 0 means fast is not used. Not implemented for MGC yet.
%%        max specifies whether to compute all 1d pairwise statistics, and the number of maximum statistics to use.
%%            0 uses the distance matrix of all dimensions, 1 uses the maximum pairwise stat, p uses the maximum p statistics averaged, etc.
%%
%% @return A list contains the following output:
%% @return corr is the sample MGC statistic within [-1,1];
%% @return localCor consists of all local correlations by double matrix index;
%% @return optimalScale the estimated optimal scale in matrix single index, which can be changed to double indices by [k,l]=ind2sub(size(localCor),optimalScale);
%% @return corrAll is all the pairwise 1D correlation when optionComponent is not 0.
%%
%% @export
%%
function [corr, localCor, optimalScale,R,corrAll]=MGC(X,Y,opts)
% check input
if nargin < 3
    opts = struct('metric','euclidean','center','mgc','fast',1,'max',0); % default parameters
end
if isfield(opts,'metric'); optionMetric = opts.metric; else optionMetric = 'euclidean'; end
if isfield(opts,'center'); optionCenter = opts.center; else optionCenter = 'mgc'; end
if isfield(opts,'fast'); optionFast = opts.fast; else optionFast = 1; end
if isfield(opts,'max'); optionMax = opts.max; else optionMax = 0; end
% check whether the given data is already kernel or distance or not
[X,indX]=checkDist(X);
[Y,indY]=checkDist(Y);
if indX>0 || indY>0
    optionFast=0;
    optionMax=0;
end
% check whether it is doing two-sample or independence testing, concatenate
% properly if it is two-sample
[X,Y]=checkTest(X,Y); 
p=size(X,2);
q=size(Y,2);

% compute test statistic
corrAll=0;
if optionMax>0
    corrAll=zeros(p,q);
    for i=1:p
        for j=1:q
            corrAll(i,j)=MGCStat(X(:,i),Y(:,j),optionMetric,optionCenter,optionFast);
        end
    end
    tmp=reshape(corrAll,p*q,1);
    corr=mean(maxk(tmp,optionMax));
    localCor=0;optimalScale=0;R=0;
else
    [corr, localCor, optimalScale,R]=MGCStat(X,Y,optionMetric,optionCenter,optionFast);
end

function [corr, localCor, optimalScale,R]=MGCStat(X,Y,optionMetric,optionCenter,optionFast)
localCor=MGCLocalCor(X,Y,optionMetric,optionCenter); % compute all localCor
[m,n]=size(localCor);
%     corr=mean(mean(localCor(2:end,2:end)));
%     %corr=localCor(m,n);
%     optimalScale=m*n;
%     R=1;
if m==1 || n==1
    corr=localCor(end);
    optimalScale=m*n;
    R=1;
else
    % find the maximal within the significant region, return optimal scale and a connectec region of local correlations R
    [corr,optimalScale, R]=MGCSmoothing(localCor,m,n); 
end

