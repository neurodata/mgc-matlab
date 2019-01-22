%% Fast Dcorr statistic that runs in O(rnlog(n)), based on:
%% C. Huang and X. Huo, "A Statistically and Numerically Efficient Independence Test based on Random Projections and Distance Covariance", 
%% https://arxiv.org/abs/1701.06054
%% A. Chaudhuri and W. Hu, “A fast algorithm for computing distance correlation,” https://arxiv.org/abs/1810.11332, 2018
%% Note that it only works for Euclidean distance.
%%
%% @param x is an n*p data matrix;
%% @param y is an n*q data matrix;
%% @param ppx is an integer specifying the number of random projection for x
%% @param ppy is an integer specifying the number of random projection for y
%% @return the mean of random projected distance correlation
%%
%% @export
%%

function [corr,ppx,ppy]=FastStatDCorRP(x,y,ppx,ppy)

if nargin<4
    ppy=100;
end
if nargin<3
    ppx=100;
end
dx=size(x,2);
dy=size(y,2);
if dx==1
    ppx=1;
end
if dy==1
    ppy=1;
end
iter=max([ppx(1),ppy(1),size(ppx,2),size(ppy,2)]);
if size(ppx,2)<2 || size(ppx,1)~=size(x,2)
    ppx=unifrnd(0,1,size(x,2),iter);
end
if dx==dy
    ppy=ppx;
else
   if size(ppy,2)<2 || size(ppy,1)~=size(y,2)
      ppy=unifrnd(0,1,size(y,2),iter);
   end
end

corr=zeros(iter,1);
for j=1:iter
    ppx(:,j)=ppx(:,j)/norm(ppx(:,j));
    ppy(:,j)=ppy(:,j)/norm(ppy(:,j));
    px=x*ppx(:,j);
    py=y*ppy(:,j);
    cov=DCovFast(px,py);
    varX=DCovFast(px,px);
    varY=DCovFast(py,py);
    corr(j)=cov/sqrt(varX*varY);
end
corr=mean(corr);
% [corr2,ind]=max(corr);
% [~,ind1]=maxk(-vecnorm(ppx-ppx(:,ind)),3)
% [~,ind2]=maxk(corr,3)
% vecnorm(ppx(:,ind2)-ppx(:,ind))
% corr3=max(corr);

