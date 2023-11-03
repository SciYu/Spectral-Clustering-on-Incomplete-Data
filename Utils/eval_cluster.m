function score = eval_cluster(label, Y)
% function score = eval_cluster(label, Y)
%
% @param label        Predicted labels using clustering method
% @param Y            True labels
% 
% @return score       Clustering scores [ACC,NMI,PUR,ARI]

baseline_missing = fieldnames(label);
for i = 1 : length(baseline_missing)
    name = baseline_missing{i};
    eval(['niter = size(label.',name,',2);']);
    for j = 1 : niter
        eval(['score.',name, '(j,:) = eval_cluster_metrics(label.', name,'(:,j), Y);']);
    end
end

end