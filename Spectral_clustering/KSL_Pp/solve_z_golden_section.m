% use golden section method
function z_ij = solve_z_golden_section(lambda, p, rho, gamma, c_ij, u_ij, alpha, max_gs_iter, gs_tol, z0)
    % Input:
    %   lambda, p, rho, gamma: hyperparameters
    %   c_ij, u_ij: variables
    %   alpha: step size
    %   max_gs_iter: maximum number of Golden Section iterations
    %   gs_tol: Golden Section method tolerance
    %   z0: initial value of z_ij
    
    % Output:
    %   z_ij: solution of z_ij
    
    a = 0;        % Lower bound of z_ij
    b = 1;        % Upper bound of z_ij
    phi = (sqrt(5) + 1) / 2;  % Golden ratio
    
    x1 = a + (1 - 1/phi) * (b - a);
    x2 = a + 1/phi * (b - a);
    
    f1 = lambda * p * x1^(p-1) + gamma * (x1 - 1) + (1/rho) * (x1 - c_ij) + u_ij;
    f2 = lambda * p * x2^(p-1) + gamma * (x2 - 1) + (1/rho) * (x2 - c_ij) + u_ij;
    
    for k = 1:max_gs_iter
        if abs(b - a) < gs_tol
            break;  % Convergence achieved
        end
        
        if f1 < f2
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (1 - 1/phi) * (b - a);
            f1 = lambda * p * x1^(p-1) + gamma * (x1 - 1) + (1/rho) * (x1 - c_ij) + u_ij;
        else
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + 1/phi * (b - a);
            f2 = lambda * p * x2^(p-1) + gamma * (x2 - 1) + (1/rho) * (x2 - c_ij) + u_ij;
        end
    end
    
    z_ij = (a + b) / 2;
end


