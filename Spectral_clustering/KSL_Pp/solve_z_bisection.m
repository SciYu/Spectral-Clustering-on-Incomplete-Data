% use bisection method to solve z
function z_ij = solve_z_bisection(lambda, p, rho, gamma, c_ij, u_ij, alpha, max_bisection_iter, bisection_tol, z0)
    % Input:
    %   lambda, p, rho, gamma: hyperparameters
    %   c_ij, u_ij: variables
    %   alpha: step size
    %   max_bisection_iter: maximum number of bisection iterations
    %   bisection_tol: bisection method tolerance
    %   z0: initial value of z_ij
    
    % Output:
    %   z_ij: solution of z_ij
    
    z_low = 0;  % Lower bound of z_ij
    z_high = 1; % Upper bound of z_ij
    
    z_ij = z0;  % Initialize z_ij
    
    for k = 1:max_bisection_iter
        f = lambda * p * z_ij^(p-1) + gamma * (z_ij - 1) + (1/rho) * (z_ij - c_ij) + u_ij;
        
        if abs(f) <= bisection_tol
            break;  % Convergence achieved
        end
        
        if f > 0
            z_high = z_ij;
        else
            z_low = z_ij;
        end
        
        z_ij = (z_low + z_high) / 2;
    end
end
