%% Compute distance or kernel correlation. Runs in O(n^2) in general. 
%% See G. Szekely and M. Rizzo and N. Bakirov, "Measuring and Testing Independence by Correlation of Distances";
%% A. Gretton and L. Gyorfi, "Consistent Nonparametric Tests of Independence";
%% C. Shen and J. Vogelstein, "The Exact Equivalence of Distance and Kernel Methods in Hypothesis Testing".
%%
%% For 1D data with Euclidean metric, utilizes a fast code that runs in O(nlogn), see
%% A. Chaudhuri and W. Hu, "A fast algorithm for computing distance correlation";
%% C. Shen and J. Vogelstein, "The Chi-Square Test of Distance Correlation".
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix for independence testing, or a m*p data matrix for two sample testing.
%% @param opts - input option structure:
%%        metric is a string that specifies which metric to use, including 'euclidean','hsic', and other variants.
%%        center is a string that specifies how the distance matrices are centered, including 'unbiased', 'biased', 'mantel' and 'mgc'.
%%        fast specifies whether fast stat computation is used or not. 1 means fast is used, 0 means fast is not used.
%%        max specifies whether to compute all 1d pairwise statistics, and the number of maximum statistics to use.
%%            0 uses the distance matrix of all dimensions, 1 uses the maximum pairwise stat, p uses the maximum p statistics averaged, etc.
%%
%% @return A list contains the following output:
%% @return corr is the distance correlation in [-1,1], while corrAll is all the pairwise 1D correlation when optionComponent is not 0.
%%
%% @export
%%
function [corr,corrAll] = DCor(X,Y,opts)
% check input
if nargin < 3
    opts = struct('metric','euclidean','center','unbiased','fast',1,'max',0); % default parameters
end
if isfield(opts,'metric'); optionMetric = opts.metric; else optionMetric = 'euclidean'; end
if isfield(opts,'center'); optionCenter = opts.center; else optionCenter = 'unbiased'; end
if isfield(opts,'fast'); optionFast = opts.fast; else optionFast = 1; end
if isfield(opts,'max'); optionMax = opts.max; else optionMax = 0; end
% check whether optionFast is eligible
if ~strcmpi(optionMetric,'euclidean') && ~strcmpi(optionMetric,'dcor')
    optionFast=0;
end
if (size(X,2)>1 || size(Y,2)>1) && optionMax==0
    optionFast=0;
end
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
if optionMax>0
    corrAll=zeros(p,q);
    for i=1:p
        for j=1:q
            corrAll(i,j)=DCorStat(X(:,i),Y(:,j),optionMetric,optionCenter,optionFast);
        end
    end
    tmp=reshape(corrAll,p*q,1);
    corr=mean(maxk(tmp,optionMax));
else
    [corr,~,~,corrAll]=DCorStat(X,Y,optionMetric,optionCenter,optionFast);
end

function [corr,varX,varY,cov]=DCorStat(X,Y,optionMetric,optionCenter,optionFast)
if optionFast==0
    % General O(n^2) version
    X=DCorInput(X,optionMetric);
    Y=DCorInput(Y,optionMetric);
    [A,B]=DCorTransform(X,Y,optionCenter,0);
    [cov]=GlobalCov(A,B'); % compute corr statistic
    [varX]=GlobalCov(A,A'); % compute varX statistic
    [varY]=GlobalCov(B,B'); % compute varY statistic
else
    % Fast O(nlogn) version for 1D Euclidean distance
    cov=DCovFast(X,Y,optionCenter);
    varX=DCovFast(X,X,optionCenter);
    varY=DCovFast(Y,Y,optionCenter);
end
corr=cov./real(sqrt(varX*varY')); % normalized
% corr=corr/size(X,1)^2;% un-normalized

if varX<=0 || varY<=0
    corr=0; % if either variance is non-positive, set corr to 0
end

% % use absolute value for mantel or Pearson correlation
% if strcmpi(optionCenter,'mantel')==true || strcmpi(optionMetric,'euclidean2')==true
%     corr=abs(corr); 
% end

function [cov]=GlobalCov(A,B)
cov=sum(sum(A.*B));
% cov=eigs(A*B*B*A,1)^0.5;