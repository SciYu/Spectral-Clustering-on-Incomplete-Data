function [Stat, Stat_std] = statistic_cluster(Score, error, error_std, recall, recall_std)
% function [stat, stat_std] = statistic_cluster(score, error, error_std, recall, recall_std)
%
% @param score              Clustering performance [ACC,NMI,PUR,ARI]
% 
% @return Stat              Statistic results [ACC,NMI,PUR,ARI,RMSE,RE,Recall]
% @return Stat_std          Statistic results with standard deviation

baseline_cluster = fieldnames(Score);
N_cluster = length(baseline_cluster);
name = baseline_cluster{1};
eval(['baseline_miss = fieldnames(Score.',name,'{1});']);
eval(['niter = length(Score.',name,');']);
N_miss = length(baseline_miss);

for i = 1 : N_cluster
    method_cluster = baseline_cluster{i};
    eval(['score = Score.',method_cluster,';']);
    stat_one = zeros(niter, N_miss, 3);
    for j = 1 : N_miss
        method_miss = baseline_miss{j};
        for k = 1 : niter
            eval(['stat_one(k,j,1) = mean(score{k}.',method_miss,'(:,1));']); % ACC
            eval(['stat_one(k,j,2) = mean(score{k}.',method_miss,'(:,2));']); % NMI
            eval(['stat_one(k,j,3) = mean(score{k}.',method_miss,'(:,3));']); % PUR
            eval(['stat_one(k,j,4) = mean(score{k}.',method_miss,'(:,4));']); % ARI
        end
    end
    eval(['stat.',method_cluster,'(1,:) = mean(stat_one(:,:,1));']); % ACC
    eval(['stat.',method_cluster,'(2,:) = mean(stat_one(:,:,2));']); % NMI
    eval(['stat.',method_cluster,'(3,:) = mean(stat_one(:,:,3));']); % PUR
    eval(['stat.',method_cluster,'(4,:) = mean(stat_one(:,:,4));']); % ARI

    eval(['stat_std.',method_cluster,'(1,:) = std(stat_one(:,:,1));']); % ACC
    eval(['stat_std.',method_cluster,'(2,:) = std(stat_one(:,:,2));']); % NMI
    eval(['stat_std.',method_cluster,'(3,:) = std(stat_one(:,:,3));']); % PUR
    eval(['stat_std.',method_cluster,'(4,:) = std(stat_one(:,:,4));']); % ARI
    
    eval(['stat_list{i} = stat.',method_cluster,';']);
    eval(['stat_std_list{i} = stat_std.',method_cluster,';']);
end

Stat = recall; Stat_std = recall_std;
rowname = {'Recall'};
for i = 1 : N_cluster
    Stat = [Stat; stat_list{i}];
    Stat_std = [Stat_std; stat_std_list{i}];
    method_cluster = baseline_cluster{i};
    rowname = [rowname, strcat(method_cluster, '-ACC'), strcat(method_cluster, '-NMI'),...
               strcat(method_cluster, '-PUR'), strcat(method_cluster, '-ARI')];
end
Stat = [Stat; error]; Stat_std = [Stat_std; error_std];
rowname = [rowname,'RE-D','RMSE-D','RE-K','RMSE-K'];

%{'Recall','ACC','NMI','PUR','ARI','RE-D','RMSE-D','RE-K','RMSE-K'});
Stat = array2table(Stat, 'VariableNames', baseline_miss, 'RowNames', rowname); 
Stat_std = array2table(Stat_std, 'VariableNames', baseline_miss, 'RowNames', rowname); 

Stat = Stat([2:end-4,end,end-2,end-1,end-3,1], [2,7:end,6,3,4,5]);
Stat_std = Stat_std([2:end-4,end,end-2,end-1,end-3,1], [2,7:end,6,3,4,5]);
end