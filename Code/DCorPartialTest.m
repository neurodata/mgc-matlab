%% Compute partial distance correlation and its p-value, given three sample data of size n*p.
%% The partial correlation and the p-value computation are based on 
%% G. Szekely and M. Rizzo, "Partial distance correlation with methods for dissimilarities";
%% C. Shen and J. Vogelstein, "The Chi-Square Test of Distance Correlation".
%%
%% @param X is an n*n distance matrix or a n*p data matrix;
%% @param Y is an n*n distance matrix or a n*q data matrix, or a m*p data matrix for two sample testing.
%% @param opts - input option structure:
%%        metric is a string that specifies which metric to use, including 'euclidean','hsic', and other variants.%% @param optionCenter is a string that specifies how the distance matrices are centered, including 'mgc'(default), 'unbiased', 'biased' and 'mantel'.
%%        max specifies whether to compute all 1d pairwise statistics, and the number of maximum statistics to use.
%%            0 uses the distance matrix of all dimensions, 1 uses the maximum pairwise stat, pq uses the maximum pq statistics averaged.
%%            Currently it only support 0,1,pq, any other integer in (1,pq) defaults to 1.
%%
%% @return A list contains the following output:
%% @return The partial correlation and its p-value
%%
%% @export

function [corr,pval] = DCorPartialTest(X,Y,Z)

n=size(X,1);
a=DCor(X,Y)-DCor(X,Z)*DCor(Z,Y);
b=(1-DCor(X,Z)^2)^0.5*(1-DCor(Y,Z)^2)^0.5;
if b==0
    corr=0;
else
    corr=a/b;
end
pval=1-chi2cdf(corr*n+1,1); 

% eps=0.000001;
% 
% [~,AC]=DCor(X,Z);
% [~,BC]=DCor(Y,Z);
% [~,CC]=DCor(Z,Z);
% 
% [A]=DCorInput(X);
% [B]=DCorInput(Y);
% [C]=DCorInput(Z);
% [A,B]=DCorTransform(A,B,'unbiased',0);
% [C,~]=DCorTransform(C,A,'unbiased',0);
% 
% PA=A-AC/CC*C;
% PB=B-BC/CC*C;
% cov=GlobalCov(PA,PB)/n/(n-3);
% varPA=GlobalCov(PA,PA)/n/(n-3);
% varPB=GlobalCov(PB,PB)/n/(n-3);
% corr=cov./real(sqrt(varPA*varPB)); % normalized
% 
% if varPA<=eps || varPB<=eps
%     corr=0; % if either variance is non-positive, set corr to 0
% end
% pval=1-chi2cdf(corr*n+1,1); 
% 
% function [cov]=GlobalCov(A,B)
% cov=sum(sum(A.*B));