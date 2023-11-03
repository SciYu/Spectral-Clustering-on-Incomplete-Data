function [imputedX] = impute_svt(Xmiss)
% function [imputedX] = impute_svt(Xmiss)
%
% Impute a data matrix. Each column is a sample. The complete matrix is obtained
% by low-rank matrix completion with the singular value thresholding (SVT) algorithm
% (see reference)
%
% @param  Xmiss       Incomplete data matrix
% @return imputedX    Imputed data matrix
%
% <Reference>
% Jian-Feng Cai, Emmanuel J Cand`es, and Zuowei Shen. "A singular value thresholding 
% algorithm for matrix completion." SIAM Journal on Optimization, 20(4):1956–1982, 2010.

low = min(Xmiss(:)); high = max(Xmiss(:));

M = Xmiss;
[n1, n2] = size(M);
r = 20;
df = r*(n1+n2-r);
m = min(5*df, round(.99*n1*n2));
p = m / (n1+n2);
I = ~isnan(Xmiss);
Omega = find(I(:));
data = M(Omega);
tau = 5*sqrt(n1*n2); 
delta = 1.2/p;    

%{
 if n1 and n2 are very different, then
   tau should probably be bigger than 5*sqrt(n1*n2)

 increase tau to increase accuracy; decrease it for speed

 if the algorithm doesn't work well, try changing tau and delta
   i.e. if it diverges, try a smaller delta (e.g. delta < 2 is a 
   safe choice, but the algorithm may be slower than necessary).
%}

maxiter = 500; 
tol = 1e-4;

[U,S,V,numiter] = SVT([n1 n2],Omega,data,tau,delta,maxiter,tol);
    
imputedX = U*S*V';
imputedX(I) = Xmiss(I);

imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;

end
