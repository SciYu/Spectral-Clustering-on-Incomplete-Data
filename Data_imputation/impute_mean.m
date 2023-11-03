function [imputedX] = impute_mean(Xmiss)
% function [imputedX] = impute_mean(Xmiss)
%
% Replace NaN values by row means.
%
% @param Xmiss      Incomplete data matrix, each column is a sample
% 
% @return imputedX  Imputed matrix with all data samples

low = min(Xmiss(:)); high = max(Xmiss(:));

idx = isnan(Xmiss);
imputedX = Xmiss;
imputedX(idx) = 0;
imputedX = imputedX + nanmean(Xmiss,2).*idx;

imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;

end