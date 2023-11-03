function [imputedX] = impute_em(Xmiss)
% function [imputedX] = impute_em(Xmiss)
%
% Impute a data matrix with the Expectation-maximization (EM) algorithm. 
%
% @param Xmiss      Incomplete data matrix, each column is a sample
% 
% @return imputedX  Imputed matrix with all data samples
%
% <Reference>
% Todd K Moon. "The expectation-maximization algorithm." 
% IEEE Signal Processing Magazine, 13(6):47–60, 1996.

low = min(Xmiss(:)); high = max(Xmiss(:));

Xmiss(isnan(Xmiss)) = -999;
imputedX = fill_with_em(Xmiss')';
        
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;

end
