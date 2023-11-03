function [K_list, D_list] = generate_k(X_list, metric, baseline_calibrate)
% function [K_list, D_list] = generate_k(X_list, metric, baseline_calibrate)
%
% @param X_list              A list of complete data matrix, each column is a sample
% @param metric              Kernel type (Gaussian, Laplacian)
% @param baseline_calibrate  Baselines for distance calibration
% 
% @return K_list             Approximate kernel matrices using baselines
% @return D_list             Approximate Euclidean distance matrices using baselines

if (nargin < 3)
    baseline_calibrate = [];
end

%% Incomplete data
Xtrue = X_list.true;
Xmiss = X_list.miss;

%% Kernel Calculation
fprintf(['\n  ', metric, ': ']);
if ~isempty(baseline_calibrate)
    %% Kernel Correction
    if strcmp(metric, 'Gaussian') || strcmp(metric, 'Laplacian')
        % Euclidean distance
        fprintf('true, ');
        D_true = distance(Xtrue, 'euclidean');
        D_miss = distance(Xmiss, 'euclidean', 'miss');
        % Distance calibration
        fprintf('dc, '); D_dc = calibrate_dc(D_miss); 
        fprintf('trf, '); D_trf = calibrate_trf(D_miss); 
        fprintf('ee, '); D_ee = calibrate_ee(D_miss);
        fprintf('kc, '); D_kc = correct_distance(D_miss); % our method
        % Euclidean distance
        D_list.true = D_true; 
        D_list.miss = D_miss;
        D_list.trf = D_trf; 
        D_list.ee = D_ee; 
        D_list.kc = D_kc;
        D_list.dc = D_dc; 
        if strcmp(metric, 'Gaussian')
            % Gaussian kernel
            K_list.true = exp(-D_true.^2 / median(D_true(:))^2);
            K_list.miss = exp(-D_miss.^2 / median(D_miss(:))^2);
            K_list.trf = exp(-D_trf.^2 / median(D_trf(:))^2);
            K_list.ee = exp(-D_ee.^2 / median(D_ee(:))^2);
            K_list.kc = exp(-D_kc.^2 / median(D_kc(:))^2);
            K_list.dc = exp(-D_dc.^2 / median(D_dc(:))^2);
        elseif strcmp(metric, 'Laplacian')
            % Laplacian kernel
            K_list.true = exp(-D_true / median(D_true(:)));
            K_list.miss = exp(-D_miss / median(D_miss(:)));
            K_list.trf = exp(-D_trf / median(D_trf(:)));
            K_list.ee = exp(-D_ee / median(D_ee(:)));
            K_list.kc = exp(-D_kc / median(D_kc(:)));
            K_list.dc = exp(-D_dc / median(D_dc(:)));
        end
    end
    
    %% Kernel Calculation
    baseline_impute = fieldnames(X_list);
    for i = 1 : length(baseline_impute)
        method = baseline_impute{i};
        if strcmp(method, 'true') || strcmp(method, 'miss')
            continue;
        else
            fprintf([method, ', ']); 
            eval(['X_impute = X_list.', method,';']);
            D_impute = distance(X_impute, 'euclidean');
            eval(['D_list.', method, ' = D_impute;']);
            switch metric
                case 'Gaussian'
                    K_impute = exp(-D_impute.^2 / median(D_impute(:))^2);
                case 'Laplacian'
                    K_impute = exp(-D_impute / median(D_impute(:)));
            end
            eval(['K_list.', method, ' = K_impute;']);
        end
    end
end

end