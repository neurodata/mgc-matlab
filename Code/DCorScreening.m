function [ind,pval,corr] = DCorScreening(X,Y,thres)
if nargin<3
    thres=0.05;
end
d=size(X,2);
pval=zeros(d,1);
corr=zeros(d,1);
ind=zeros(d,1);
for i=1:d
    [tmp,corr(i)]=DCorFastTest(X(:,i),Y);
    pval(i)=tmp;
    if tmp<thres
        ind(i)=1;
    end
end