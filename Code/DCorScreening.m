function [ind,corr] = DCorScreening(X,Y,option)

% if nargin<3
%     thres=0.05;
% end
% d=size(X,2);
% pval=zeros(d,1);
% corr=zeros(d,1);
% for i=1:d
%     [corr(i),pval(i)]=DCorFastTest(X(:,i),Y);
% end
% ind=(pval<thres);

if nargin<3
    option=0;
end

if option==0
    opts = struct('center','unbiased');
else
    opts = struct('center','biased');
end

[n,d]=size(X);
thres=5/n^0.9;
corr=zeros(d,1);
for i=1:d
    if option==2
        tmp=corrcoef(X(:,i),Y);
        corr(i)=abs(tmp(1,2));
    else
        corr(i)=DCor(X(:,i),Y,opts);
    end
end

ind=(corr>thres);