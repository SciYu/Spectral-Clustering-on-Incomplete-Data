function [imputedX] = impute_gr(Xmiss)
% function [imputedX] = impute_gr(Xmiss)
%
% Impute a data matrix. Each column is a sample. Each NaN value is replaced
% by the statistical value calculated by the GR method (see reference).
%
% @param Xmiss      Incomplete data matrix, each column is a sample
% @return imputedX  Imputed matrix with all data samples
%
% <Reference>
% Balzano, Laura, Robert Nowak, and Benjamin Recht. "Online identification 
% and tracking of subspaces from highly incomplete information." 2010 48th 
% Annual allerton conference on communication, control, and computing (Allerton). IEEE, 2010.

low = min(Xmiss(:)); high = max(Xmiss(:));

max_rank = 5; 
step_size = 0.1;
max_Cycles = 50;

Xzero = impute_zero(Xmiss);
M = double(~isnan(Xmiss));
[numr, numc] = size(Xmiss);
[I, J, S] = miss_grouse(Xmiss, M);
[U,V,~] = grouse(I,J,S,numr,numc,max_rank,step_size,max_Cycles);
imputedX = U*V';
imputedX = Xzero + imputedX .* (1-M);

imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
imputedX(isnan(imputedX)) = 0;

end

%%
function [I, J, S] = miss_grouse(Xmiss, M)

I = []; J = []; S = []; kk = 0;
[numr, numc] = size(Xmiss);
for ii = 1:numc
    for jj = 1:numr
        if M(jj,ii) == 1
            kk = kk + 1;
            J(kk,1) = ii;
            I(kk,1) = jj;
            S(kk,1) = Xmiss(jj,ii);
        end
    end
end
end