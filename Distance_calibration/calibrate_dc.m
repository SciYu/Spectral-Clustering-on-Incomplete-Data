function [D_new] = calibrate_dc(D0)
% function [D_new] = calibrate_dc(D0)
%
% Calibrate a distance matrix by the Double-Centering (DC) algorithm.
% 
% @param D0      Initial distance matrix (estimated from incomplete data)
%
% @return D_new  Calibrated distance matrix
%
% <Reference>
% Gisbrecht, A.; and Schleif, F.-M. 2015. Metric and non-metric proximity 
% transformations at linear costs. Neurocomputing, 167: 643–657.

%% Step 1. Transfer distance matrix to similariy matrix
n = size(D0,1);
J = eye(n) - (1/n) * ones(n); % centering matrix

% Apply double centering
S = -0.5 * J * (D0.^2) * J;

%% Step 2. Obtain a PSD similarity matrix
[U, Lambda, ~] = svd(S);

% Clip-operation: set all negative eigenvalues to zero
Lambda_new = max(real(Lambda), 0);

% Obtain a PSD similarity matrix
S_new = real(U * Lambda_new * U');
S_new = (S_new + S_new') / 2;

%% Step 3. Transfer similarity matrix to distance matrix
D_new = transfer_StoD(S_new);

end


%%
function [D] = transfer_StoD(S)
% D_ij^2 = S_ii + S_jj - 2*S_ij
n = size(S,1);
s = diag(S);
S_ii = repmat(s,1,n);
S_jj = repmat(s',n,1);
D = sqrt(S_ii + S_jj - 2*S);
end
