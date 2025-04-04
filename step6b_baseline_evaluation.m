%% This code can be run directly after changing file
clear
file = 'E:\Zhongwei\SAH Code Publish\models\Baseline_RF\rf_50_trial10_seed19284_allFea.mat';
rf_para = 50;
trial_num = 10;

load(file);

%% calculate mean AUC of all trials

time_step = 6;
time_range = 48:time_step:240;
num_fold = 5;
confidence_interval = [];
for itrials = 1:length(preds_all)
    labels = labels_all{itrials};
    pred_probability = pred_probability_all{itrials};
    preds = preds_all{itrials};
    for i = 1:length(time_range)
        label_tmp = [];
        pred_prob_tmp = [];
        preds_tmp = [];
        acc_tmp = [];
        auc_tmp = [];
        for j = 1:num_fold
            label_tmp = [label_tmp;labels{i,j}];
            pred_prob_tmp = [pred_prob_tmp;pred_probability{i,j}];
            preds_tmp = [preds_tmp;preds{i,j}];
        end
        acc_all(itrials,i) = sum(label_tmp==cell2mat(preds_tmp))/length(preds_tmp);
        [X_tmp,Y_tmp,T_tmp,AUC_tmp] = perfcurve(label_tmp,1-pred_prob_tmp(:,1),true);
        AUC_all(itrials,i) = AUC_tmp;
        [A, Aci] = auc([label_tmp,1-pred_prob_tmp(:,1)],0.05,'hanley');
        confidence_interval{itrials,i} = Aci;
    end
end
acc = acc_all;
acc_std = zeros(size(acc));
AUC = AUC_all;
AUC_std = zeros(size(AUC));
mean_trial_AUC = mean(AUC);

%%
figure
hold on,
plot(1:length(time_range),mean_trial_AUC,'-b'); %,'linewidth',2
set(gca, 'XTick', 1:4:length(time_range), 'XTickLabel', strsplit(num2str(2:10)));
xlabel('Days After Bleed');
ylabel('AUC');
title('Mean AUC of DCI Prediction')
xlim([0 length(time_range)+1])
%%
figure,
x = 1:length(time_range);
fr = fit(x',mean_trial_AUC','smoothingspline','SmoothingParam',0.1);
plot(fr,'-r');
set(gca, 'XTick', 1:4:length(time_range), 'XTickLabel', strsplit(num2str(2:10)));
xlabel('Days After Bleed');
ylabel('AUC');
title('Mean AUC of DCI Prediction')
xlim([0 length(time_range)+1])
hold on,scatter(x,AUC,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
legend('selected+cumulative')



%% calculate mean AUC - average all past results of all trials
% current time range/model result will take past prediction into consideration

time_step = 6;
time_range = 48:time_step:240;
num_fold = 5;
confidence_interval = [];
for itrials = 1:length(preds_all)
    labels = labels_all{itrials};
    pred_probability = pred_probability_all{itrials};
    preds = preds_all{itrials};
    label_tmp = [];
    pred_prob_tmp = [];
    preds_tmp = [];
    for i = 1:length(time_range)
        
        % acc_tmp = [];
        auc_tmp = [];
        for j = 1:num_fold
            label_tmp = [label_tmp;labels{i,j}];
            pred_prob_tmp = [pred_prob_tmp;pred_probability{i,j}];
            preds_tmp = [preds_tmp;preds{i,j}];
        end
        % acc_all(itrials,i) = sum(label_tmp==cell2mat(preds_tmp))/length(preds_tmp);
        [X_tmp,Y_tmp,T_tmp,AUC_tmp] = perfcurve(label_tmp,1-pred_prob_tmp(:,1),true);
        AUC_all(itrials,i) = AUC_tmp;
        [A, Aci] = auc([label_tmp,1-pred_prob_tmp(:,1)],0.05,'hanley');
        confidence_interval{itrials,i} = Aci;
    end
end
% acc = acc_all;
% acc_std = zeros(size(acc));
AUC = AUC_all;
AUC_std = zeros(size(AUC));
mean_trial_AUC = mean(AUC);

figure
hold on,
plot(1:length(time_range),mean_trial_AUC,'-b'); %,'linewidth',2
set(gca, 'XTick', 1:4:length(time_range), 'XTickLabel', strsplit(num2str(2:10)));
xlabel('Days After Bleed');
ylabel('AUC');
title('Mean AUC of DCI Prediction - Past Considered')
xlim([0 length(time_range)+1])

figure,
x = 1:length(time_range);
fr = fit(x',mean_trial_AUC','smoothingspline','SmoothingParam',0.1);
plot(fr,'-r');
set(gca, 'XTick', 1:4:length(time_range), 'XTickLabel', strsplit(num2str(2:10)));
xlabel('Days After Bleed');
ylabel('AUC');
title('Mean AUC of DCI Prediction - Past Considered')
xlim([0 length(time_range)+1])
hold on,scatter(x,AUC,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
legend('selected+cumulative')


%% average the past time range AUC

indices = 1:length(mean_trial_AUC);
average_auc = cumsum(mean_trial_AUC) ./ indices;
figure
hold on,
plot(1:length(time_range),average_auc,'-b'); %,'linewidth',2
set(gca, 'XTick', 1:4:length(time_range), 'XTickLabel', strsplit(num2str(2:10)));
xlabel('Days After Bleed');
ylabel('AUC');
title('Mean AUC of DCI Prediction - Past Average')
xlim([0 length(time_range)+1])
%%
figure,
x = 1:length(time_range);
fr = fit(x',average_auc','smoothingspline','SmoothingParam',0.1);
plot(fr,'-r');
set(gca, 'XTick', 1:4:length(time_range), 'XTickLabel', strsplit(num2str(2:10)));
xlabel('Days After Bleed');
ylabel('AUC');
title('Mean AUC of DCI Prediction - Past Average')
xlim([0 length(time_range)+1])
hold on,scatter(x,AUC,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
legend('selected+cumulative')
