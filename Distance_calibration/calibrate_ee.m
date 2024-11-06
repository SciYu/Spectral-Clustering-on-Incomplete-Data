function [D_new] = calibrate_ee(D0, epsilon, maxiter)
% function [D_new] = calibrate_ee(D0, epsilon, maxiter)
%
% Calibrate a distance matrix by a Euclidean Embedding (EE) method. (see reference)
%
% @param  D0        Initial distance matrix (estimated from incomplete data)
% @param  epsilon   Default 0.02
% @param  maxiter   Maximum iterations (default 100)
%
% @return D_new     Calibrated distance matrix
%
% <Reference>
% [1] Wenye Li, Fangchen Yu. "Calibrating Distance Metrics Under Uncertainty." ECML, 2022.
% [2] Fangchen Yu, et al. "Highly-Efficient Robinson-Foulds Distance Estimation with Matrix Correction" ECAI, 2023.

if (nargin < 3)
    maxiter = 100;
end
if (nargin < 2)
    epsilon = 0.02;
end

n = size(D0,1);

gamma = -epsilon / max(D0(:));
low = exp(-epsilon); high = 1;
K0 = exp(gamma.*D0);
K = nearpsd(K0, maxiter, low, high);
D_new = log(K) ./ gamma;

D_new = (D_new+D_new') / 2;
D_new(1:n+1:n^2) = 0;
D_new = max(real(D_new), 0);

end


%%
function Y = nearpsd(A, maxits, low, high)
% function Y = nearpsd(A, maxits, low, high)
%
% Computes the nearest positive semi-definite matrix 
% for a given square matrix.
%
% @param A        a squared matrix to be calibrated
% @param maxits   max num of iters allowed, default 100
% @param low      default 0 
% @param high     default 1
%
% @return Y       nearest psd matrix to A

if  ~isequal(A,A')
    A = (A + A') / 2;
end
if nargin < 4
    high = 1;
end
if nargin < 3
    low = 0; 
end
if nargin < 2
    maxits = 100;
end

% threshold for convergence & eigs
tolconv = 1.0e-6;
toleigs = 1.0e-5;

n = size(A,1);

U = zeros(n);
Y = A;

[V, D] = eig(Y);
d = diag(D);

iter = 0;
while 1
    T = Y - U;

    % project onto psd matrices
    [Q, D] = eig(T);
    d = diag(D);
    p = d > toleigs*d(n);
    X = Q(:,p) * D(p,p) * Q(:,p)';

    % update correction
    U = X - T;

    % maximum iteration & convergence test
    iter = iter + 1;
    if iter == maxits
        %fprintf('Max iterations reached. ');
        break; 
    end
    if norm(Y-X,'inf')/norm(Y,'inf') <= tolconv 
    	break;
    end
    
    % problem-dependent here
    Y = X;
    Y(1:n+1:n*n) = 1;
    Y(Y<low) = low;
    Y(Y>high) = high;
end

Y(1:n+1:n*n) = 1;
Y(Y<low) = low;
Y(Y>high) = high;

end
