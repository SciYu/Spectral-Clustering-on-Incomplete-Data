function [imputedX] = impute_zero(Xmiss)
% function [imputedX] = impute_zero(Xmiss)
%
% Impute a data matrix. Each NaN value in a sample is replaced by zero.
%
% @param  Xmiss       Incomplete data matrix
% @return imputedX    Imputed data matrix

imputedX = Xmiss;
imputedX(isnan(Xmiss)) = 0;

end