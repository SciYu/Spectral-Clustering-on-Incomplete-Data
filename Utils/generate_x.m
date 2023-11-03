function X_list = generate_x(X, idx, baseline_impute)
% function X_list = generate_x(X, idx, baseline_impute)
%
% @param X                  Complete data matrix, each column is a sample
% @param idx                Index matrix for missing values
% @param baseline_impute    Data Imputation Baselines
% 
% @return X_list            Imputed data matrices using baselines

%% Incomplete Data
Xtrue = full(X);
Xmiss = Xtrue;
Xmiss(idx) = NaN;
X_list.true = Xtrue;
X_list.miss = Xmiss;

%% Data Imputation
for i = 1 : length(baseline_impute)
    method = baseline_impute{i};
    fprintf([method,', ']);
    eval(['X_impute = impute_', method,'(Xmiss);']);   
    X_impute(isnan(X_impute)) = 0;
    eval(['X_list.', method, ' = X_impute;']);
end

end