clear
close all
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_vas.mat');
vascular_features = vas_features_dci;

load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_asym.mat');
asymmetry = asym_features_dci;
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_global.mat');

spike = spike_dci;
features = features_dci;

%%
for rf_para = 50%20:10:100
for itrials = 1%:20
time_step = 6;
time_range = (-5*24:time_step:3*24)+1000;
num_fold = 5;

preds = cell(length(time_range),num_fold);
labels = cell(length(time_range),num_fold);
pred_probability = cell(length(time_range),num_fold);
weights = cell(length(time_range),num_fold);
for i = 1:length(time_range)
    pts_indicators = [];
    X_all = [];
    Y_all = [];
    
    for iblock = 1:time_step
        x = [];
        spike_x = spike(:,time_range(i)-iblock+1)';
        x = [spike_x]; 

        %% selected features
        for ifea = [2,7] %% adr pav
            feature_tmp = features{ifea}(:,time_range(i)-iblock+1)';
            x = [x; feature_tmp];
        end
        for ifea = [2,7] %% adr pav
            for jfea = 1:size(vascular_features,2)
                feature_tmp = vascular_features{ifea,jfea}(:,time_range(i)-iblock+1)';
                if ifea==3||ifea==4||ifea==5||ifea==6
                    feature_tmp = 10*log10(feature_tmp);
                end
                x = [x; feature_tmp];
            end
        end
        for ifea = 6 %% alpha pca left and right
            for jfea = [3 6]
                feature_tmp = vascular_features{ifea,jfea}(:,time_range(i)-iblock+1)';
                if ifea==3||ifea==4||ifea==5||ifea==6
                    feature_tmp = 10*log10(feature_tmp);
                end
                x = [x; feature_tmp];
            end
        end
        for ifea = [3 4] %% total delta aca left and right
            for jfea = [1 4]
                feature_tmp = vascular_features{ifea,jfea}(:,time_range(i)-iblock+1)';
                if ifea==3||ifea==4||ifea==5||ifea==6
                    feature_tmp = 10*log10(feature_tmp);
                end
                x = [x; feature_tmp];
            end
        end
        for ifea = 6 % asymmetry alpha mca pca
            for jfea = [2,3]
                feature_tmp = asymmetry{ifea,jfea}(:,time_range(i)-iblock+1)';
                x = [x; feature_tmp];
            end
        end
        for ifea = 4 % asymmetry delta mca
            for jfea = [2]
                feature_tmp = asymmetry{ifea,jfea}(:,time_range(i)-iblock+1)';
                x = [x; feature_tmp];
            end
        end
        for ifea = 1 % asymmetry shannon aca mca
            for jfea = [1 2]
                feature_tmp = asymmetry{ifea,jfea}(:,time_range(i)-iblock+1)';
                x = [x; feature_tmp];
            end
        end
        for ifea = [3 5] % asymmetry theta total mca
            for jfea = [2]
                feature_tmp = asymmetry{ifea,jfea}(:,time_range(i)-iblock+1)';
                x = [x; feature_tmp];
            end
        end

        %% cumulative features adr pav spikes
%         feature_tmp = cumsum(spike(:,time_range(1)-time_step+1:time_range(i)-time_step+1),2);
%         x = [x; feature_tmp(:,end)'];
%         feature_tmp = cumsum(features{2}(:,time_range(1)-time_step+1:time_range(i)-time_step+1),2);
%         x = [x; feature_tmp(:,end)'];
%         feature_tmp = cumsum(features{7}(:,time_range(1)-time_step+1:time_range(i)-time_step+1),2);
%         x = [x; feature_tmp(:,end)'];
        
        
        index = find(~isnan(mean(x)));
        pts_indicators = [pts_indicators,index];
        X = x(:,index)';
        Y = dci_labels_new(index);

        
        X_all = [X_all;X];
        Y_all = [Y_all;Y'];
    end
        
        
    pts_unique = unique(pts_indicators);
    num_pts = length(pts_unique);
    num_test = round(num_pts/num_fold);
    idx = randperm(num_pts);
    pts_unique = pts_unique(idx);

    X_all = normalize(X_all);

    for ifold = 1:num_fold
        if ifold~=num_fold
            start_index = (ifold-1)*num_test+1;
            end_index = ifold*num_test;
        else
            start_index = (ifold-1)*num_test+1;
            end_index = num_pts;
        end
        
        idx = arrayfun( @(x)( find(pts_indicators==x) ), pts_unique(start_index:end_index),'UniformOutput',false);
        test_index = cell2mat(idx);
        train_index = setdiff(1:num_pts,start_index:end_index);
        idx = arrayfun( @(x)( find(pts_indicators==x) ), pts_unique(train_index),'UniformOutput',false);
        train_index = cell2mat(idx);
        test_data = X_all(test_index,:);
        test_label = Y_all(test_index);
        train_data = X_all(train_index,:);
        train_label = Y_all(train_index);
        
        %% Random Forest
        B=TreeBagger(rf_para,train_data,train_label,'Method','Classification','OOBVarImp','On');
        weight_tmp=B.OOBPermutedVarDeltaError;
        [pred, probabilities] = predict(B, test_data);
        pred = str2double(pred);
  
        preds{i,ifold} = pred;
        labels{i,ifold} = test_label;
        pred_probability{i,ifold} = probabilities;

    end
end
preds_all{itrials} = preds;
labels_all{itrials} = labels;
pred_probability_all{itrials} = pred_probability;


end

save(['E:\Zhongwei\SAH Code Publish\models\Baseline_RF\rf_',num2str(rf_para),'_trial',num2str(num_itrials), ...
        '_seed',num2str(seed),'_slectFea_alignDCI'],'preds_all','labels_all','pred_probability_all','time_range');
end
%% for demo
num_itrials = 1;
seed = 1;
save(['E:\Zhongwei\SAH Code Publish\models\Baseline_RF\rf_',num2str(rf_para),'_trial',num2str(num_itrials), ...
        '_seed',num2str(seed),'_slectFea_alignDCI'],'preds_all','labels_all','pred_probability_all','time_range');