% An auxiliary function to check size of X and Y, and concatenate data if
% for two-sample testing.
function [X,Y]=checkTest(X,Y)

m=size(X,1);
n=size(Y,1);
if m ~=n && size(X,2)==size(Y,2)
    %disp('The dependence measure is computed for two-sample testing')
    U=[X;Y];
    Y=[zeros(m,1);ones(n,1)];
    X=U;
else
    return
end