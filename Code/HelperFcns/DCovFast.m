%% Fast Distance Covariance that runs in O(nlog(n)) for 1D data in Euclidean distance, based on:
%% A. Chaudhuri and W. Hu, “A fast algorithm for computing distance correlation";
%% C. Shen and C. Priebe and J. Vogelstein, "From Distance Correlation to Multiscale Graph Correlation";
%% C. Shen and J. Vogelstein, "The Chi-Square Test of Distance Correlation".
%%
%% @param x is an n*1 data matrix;
%% @param y is an n*1 data matrix;
%% @param option is a string that specifies the centering, including 'mgc'(default), 'unbiased', 'biased' and 'simple'.
%% @return the distance covariance (note that this is an un-normalized covariance, see DCor.m for normalized correlation)
%%
%% @export
%%

function covsq=DCovFast(x,y,option)
if nargin < 3
    option='mgc'; % default centering scheme
end
n=length(x);

% sort the data
[x,Index]=sort(x);
y=y(Index);

% compute term1-3
si=cumsum(x);
s=si(n);
ax=(-(n-2):2:n).'.*x+(s-2*si);

v=[x y x.*y];
nw=size(v,2);

idx=zeros(n,2);
idx(:,1)=1:n;

iv1=zeros(n,1);iv2=zeros(n,1);iv3=zeros(n,1);iv4=zeros(n,1);

i=1;r=1;s=2;
while i<n
    gap=2*i;
    k=0;
    idxr=idx(:,r);
    csumv=[zeros(1,nw);cumsum(v(idxr,:))];
    for j=1:gap:n
        st1=j;e1=min(st1+i-1,n);
        st2=j+i;e2=min(st2+i-1,n);
        while (st1<=e1)&&(st2<=e2)
            k=k+1;
            idx1=idxr(st1);
            idx2=idxr(st2);
            if y(idx1)>=y(idx2)
                idx(k,s)=idx1;
                st1=st1+1;
            else
                idx(k,s)=idx2;
                st2=st2+1;
                iv1(idx2,1)=iv1(idx2)+e1-st1+1;
                iv2(idx2)=iv2(idx2)+(csumv(e1+1,1)-csumv(st1,1));
                iv3(idx2)=iv3(idx2)+(csumv(e1+1,2)-csumv(st1,2));
                iv4(idx2)=iv4(idx2)+(csumv(e1+1,3)-csumv(st1,3));
            end
        end
        if st1<=e1
            kf=k+e1-st1+1;
            idx((k+1):kf,s)=idxr(st1:e1,:);
            k=kf;
        else
            if st2<=e2
                kf=k+e2-st2+1;
                idx((k+1):kf,s)=idxr(st2:e2,:);
                k=kf;
            end
        end
    end
    i=gap;
    r=3-r;s=3-s;
end

covterm=n*(x-mean(x)).'*(y-mean(y));

c1=iv1.'*v(:,3);
c2=sum(iv4);
c3=iv2.'*y;
c4=iv3.'*x;
d=4*((c1+c2)-(c3+c4))-2*covterm;

ySorted=y(idx(n:-1:1,r));
si=cumsum(ySorted);
s=si(n);
by=zeros(n,1);
by(idx(n:-1:1,r))=(-(n-2):2:n).'.*ySorted+(s-2*si);

term1=d;
term2=(ax.'*by);
term3=sum(ax)*sum(by);

% assign weights based on different centering scheme
switch option
    case 'mgc' % almost unbiased transform utilized by mgc
        n1=n*(n-1);
        n2=n1*(n-1);
        n3=n2*(n-1);
        term3=term3-term2;
    case 'unbiased' % unbiased dcor transform
        n1=n*(n-3);
        n2=n1*(n-2);
        n3=n2*(n-1);
    case 'biased' % biased dcor transform
        n1=n*n;
        n2=n1*n;
        n3=n2*n;
    case 'mantel' % almost unbiased transform utilized by mgc
        n1=n*n;
        n2=1;
        n3=-n1*(n-1)^2*n/(n-2);
        term2=0;
    otherwise % default to mgc transform
        n1=n*(n-1);
        n2=n1*(n-1);
        n3=n2*(n-1);
        term3=term3-term2;
end
covsq=term1/n1-2*term2/n2+term3/n3;