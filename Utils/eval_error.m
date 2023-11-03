function [error, error_std] = eval_error(D_list, K_list)
% function [error, error_std] = eval_error(W_list, K_list)
%
% @param W_list             List of similarity/distance matrices
% @param K_list             List of kernel matrices
% 
% @return error             Statistic results [re, rmse]
% @return error_std         Statistic results with standard deviation

baseline = fieldnames(D_list{1});
N_miss = length(baseline);
niter = length(D_list);

for i = 1 : niter
    D = D_list{i};
    Dtrue = D.true;
    Dmiss = D.miss;
    K = K_list{i};
    Ktrue = K.true;
    Kmiss = K.miss;
    for j = 1 : N_miss
        method = baseline{j};
        eval(['D0 = D.',method,';']);
        eval(['K0 = K.',method,';']);
        
        % re = relative error
        D_re(i,j) = norm(D0-Dtrue,'fro') / norm(Dtrue,'fro');
        K_re(i,j) = norm(K0-Ktrue,'fro') / norm(Ktrue,'fro');
        
        % rmse = relative-mean-square error
        D_rmse(i,j) = norm(D0-Dtrue,'fro')^2 / norm(Dmiss-Dtrue,'fro')^2;
        K_rmse(i,j) = norm(K0-Ktrue,'fro')^2 / norm(Kmiss-Ktrue,'fro')^2;
    end
end

error = [mean(D_re); mean(D_rmse); mean(K_re); mean(K_rmse)];
error_std = [std(D_re); std(D_rmse); std(K_re); std(K_rmse)];

end