function [Recall] = eval_recall_single(W0, Wtrue, metric, topk)
% function [Recall] = eval_recall_single(W0, Wtrue, metric, topk)
%
% Evaluation function of recall performance.
%
% @param W0          Estimated weighted distance/similarity matrix
% @param Wtrue       True weighted distance/similarity matrix
% @param metric      'distance' or 'similarity'
% @param topk        Top-k similar items to calculate the recall
%
% @return recall     Recall = overlapping ratio of two sets of similar items

n = size(Wtrue, 1);
recall = zeros(1, n); 
if strcmp(metric, 'distance')
    order = 'ascend';
    Wtrue(1:n+1:n^2) = inf;
    W0(1:n+1:n^2) = inf;
elseif strcmp(metric, 'similarity')
    order = 'descend';
    Wtrue(1:n+1:n^2) = 0;
    W0(1:n+1:n^2) = 0;
end

[~, Imiss_list] = sort(W0, order);
[~, Itrue_list] = sort(Wtrue, order);
Imiss_top = Imiss_list(1:topk,:);
Itrue_top = Itrue_list(1:topk,:);
for i = 1 : n
    recall(1,i) = length(intersect(Itrue_top(:,i), Imiss_top(:,i)));
end
Recall = mean(recall / topk);

end
