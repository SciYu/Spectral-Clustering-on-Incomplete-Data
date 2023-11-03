function [D_new] = calibrate_trf(D0, tol)
% function [D_new] = calibrate_trf(D0, tol)
%
% Calibrate a distance matrix by the Triangle Fixing (TRF) algorithm (see reference)
%
% @param  D0      Initial distance matrix
% @param  tol     Tolerance for terminal
% @return D_new   Calibrated distance matrix
%
% <Reference>
% Brickell, Justin, et al. "The metric nearness problem." SIAM Journal on 
% Matrix Analysis and Applications 30.1 (2008): 375-396.

if (nargin < 2)
    tol = 1e-4;
end

D0 = triu(D0, 1);      % upper triangular part of the matrix
n = size(D0, 1);      % size of D matrix
z = zeros(n, n, n);  % dual varaiables

E = zeros(n);        % error matrix
delta = 1 + tol;     % parameter for convergence
Niter = 100;         % number of iterations
niter = 0;

while delta > tol && niter <= Niter
    niter = niter + 1;
    E0 = E;
    for i = 1 : n-2
        for k = i+1 : n-1
            for j = k+1 : n
                v = D0(i,k) + D0(k,j) - D0(i,j);
                mu = 1/3 * (E(i,j) - E(i,k) - E(k,j) - v);
                theta = max(mu, -z(i,k,j));
                E(i,j) = E(i,j) - theta;
                E(i,k) = E(i,k) + theta;
                E(k,j) = E(k,j) + theta;
                z(i,k,j) = z(i,k,j) + theta;
            end
        end
    end
    diff = abs(E - E0);
    delta = sum(diff(:));
end
D_new = D0 + E;
D_new = D_new + D_new';
D_new(1:n+1:n*n) = 0;
D_new = max(D_new, 0);

end
