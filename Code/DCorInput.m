%% Process the input into proper distance matrix or kernel induced distance matrix, depending on choice of the metric / kernel.
%%
%% @param X can be n*p data, n*n dissimilarity, or n*n similarity.
%% @param optionMetric is a string that specifies which global correlation to compute, including 'dcor'(default),'hsic', and other variants.
%%
%% @return X as the n by n processed distance matrix.
%%
%% @export
%%

function [X]=DCorInput(X,optionMetric)
if nargin<2
    optionMetric='euclidean';
end
if size(X,1)==size(X,2) && sum(diag(X))==size(X,1) % similarity matrix
    X=1-X/max(max(X));
    return;
end

if size(X,1)==size(X,2) && sum(diag(X))==0
    % dissimilarity matrix
else
    X=squareform(pdist(X));
end

% for the kernel choice, they are transformed to kernel induced distance
% matrix via 1-X/max(max(X)). As all kernel choice below are translation
% invariant and equals 1 in maximum, it is just 1-X
if strcmpi(optionMetric,'gaussian01')
    deg = max(max(X))*0.1;
    X=exp(-X.^2/deg^2);
    X=1-X;
end
if strcmpi(optionMetric,'laplace01')
    deg = max(max(X))*0.1;
    X=exp(-X/deg);
    X=1-X;
end
if strcmpi(optionMetric,'gaussianMax')
    deg = max(max(X));
    X=exp(-X.^2/deg^2);
    X=1-X;
end
if strcmpi(optionMetric,'laplaceMax')
    deg = max(max(X));
    X=exp(-X/deg);
    X=1-X;
end
if strcmpi(optionMetric,'laplace')
    deg = median(X(X>0));
    X=exp(-X/deg);
    X=1-X;
end
if strcmpi(optionMetric,'hsic')
    deg = sqrt(0.5*median(X(X>0)).^2);
    X=exp(-X.^2/2/deg^2);
    X=1-X;
end
if strcmpi(optionMetric,'euclidean01')
    X=X.^0.1;
end
if strcmpi(optionMetric,'euclidean19')
    X=X.^1.9;
end
if strcmpi(optionMetric,'euclidean2') || strcmpi(optionMetric,'pearson')
    X=X.^2;
end

% scale=0;
% if scale==1
%     X=(X-min(min(X)))/(max(max(X))-min(min(X)));
% end
