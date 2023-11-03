function [Recall, Recall_std] = eval_recall(W_list, metric, topk)
% function [Recall, Recall_std] = eval_recall(W_list, metric, topk)
%
% Evaluation function of recall performance.
%
% @param W_list       A list of weighted distance/similarity matrix
% @param metric       'distance' or 'similarity'
% @param topk         Top-k similar items to calculate the recall
%
% @return Recall      Recall = overlapping ratio of two sets of similar items
% @return Recall_std

baseline = fieldnames(W_list{1});
niter = length(W_list);
N_miss = length(baseline);

recall = zeros(niter, N_miss);
for k = 1 : niter
    Wtrue = W_list{k}.true;
    for j = 1 : N_miss
        method = baseline{j};
        eval(['W0 = W_list{k}.',method,';']);
        recall(k,j) = eval_recall_single(W0, Wtrue, metric, topk);
    end
end
Recall = mean(recall);
Recall_std = std(recall);

end
