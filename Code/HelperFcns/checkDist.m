% An auxiliary function to check whether X is a distance or kernel matrix,
% return the induced distance matrix X, and ind=0 means not distance nor
% kernel, ind=1 means distance, ind=2 means kernel
function [X,ind]=checkDist(X)

ind=0;
if size(X,1)==size(X,2) && sum(diag(X))==size(X,1) % similarity matrix
    X=1-X/max(max(X));
    ind=2;
    return;
end

if size(X,1)==size(X,2) && sum(diag(X))==0
    % dissimilarity matrix
    ind=1;
    return;
end