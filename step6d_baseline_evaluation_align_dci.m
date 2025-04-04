clear
load('E:\Zhongwei\SAH Code Publish\models\Baseline_RF\rf_50_trial1_seed1_slectFea_alignDCI.mat');
time_step = 6;
time_range = -5*24:time_step:3*24;
num_fold = 5;
%%
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
        if sum(pred_probability{i,j}(:,1)==1)==length(pred_probability{i,j}(:,1))
            pred_probability{i,j}(:,2) = 1-pred_probability{i,j}(:,1);
        end
        pred_prob_tmp = [pred_prob_tmp;pred_probability{i,j}];
        preds_tmp = [preds_tmp;preds{i,j}];

    end

    acc_all(itrials,i) = sum(label_tmp==preds_tmp)/length(preds_tmp);
    [X_tmp,Y_tmp,T_tmp,AUC_tmp] = perfcurve(label_tmp,1-pred_prob_tmp(:,1),true);
    AUC_all(itrials,i) = AUC_tmp;

end
end

acc = acc_all;
acc_std = zeros(size(acc));
AUC = AUC_all;
AUC_std = zeros(size(AUC));
max(AUC)

%%
figure,
x = 1:length(time_range);
fr = fit(x',AUC','smoothingspline','SmoothingParam',0.1);
plot(fr,'-r');
set(gca, 'XTick', 1:4:length(time_range), 'XTickLabel', strsplit(num2str(-5:3)));
xlabel('Days After DCI');
ylabel('AUC');
title('Mean AUC of DCI Prediction')
xlim([0 length(time_range)+1])
set(gca,'tickdir','out')
hold on,scatter(x,AUC,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
