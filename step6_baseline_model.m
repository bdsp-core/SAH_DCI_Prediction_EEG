clear
close all
% data used is not interpolated
% all feature is used
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_remove_dci_vas.mat')
vascular_features = features;
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_remove_dci_asym.mat')
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_remove_dci_global.mat')
%%
% for rf_para = 50%20:10:100
% for itrials = 1:20
num_itrials=1;%5;
time_step = 6;% 12;?
time_range = 48:time_step:240; % 48:time_step:288;
num_fold = 5;
seed = 19284;
rng(seed);

%%
X_all_range = cell(length(time_range),1);
Y_all_range = cell(length(time_range),1);
pts_indicators_range = cell(length(time_range),1);
for i = 1:length(time_range)
    pts_indicators = [];
    X_all = [];
    Y_all = [];
    
    % data setup: 1(spike)+7(spectral)+7*6(vascular)+7*3(asymm)=71
    for iblock = 1:time_step
        x = [];
        % spike
        spike_x = spike(:,time_range(i)-iblock+1)';
        x = [spike_x]; 
        
        % global features
        for ifea = 1:length(features)
            feature_tmp = features{ifea}(:,time_range(i)+iblock-time_step)';
            if ifea==3||ifea==4||ifea==5||ifea==6
                % 'TotalPow','DeltaPow','ThetaPow','AlphaPow'
                feature_tmp = 10*log10(feature_tmp);
            end
            x = [x; feature_tmp];
        end

        % vascular features
        for ifea = 1:size(vascular_features,1) 
            for jfea = 1:size(vascular_features,2)
                feature_tmp = vascular_features{ifea,jfea}(:,time_range(i)-iblock+1)';
                if ifea==3||ifea==4||ifea==5||ifea==6
                    feature_tmp = 10*log10(feature_tmp);
                end
                x = [x; feature_tmp];
            end
        end
        % asymmetry features
        for ifea = 1:size(asym_features,1) 
            for jfea = 1:size(asym_features,2)
                feature_tmp = asym_features{ifea,jfea}(:,time_range(i)-iblock+1)';
                x = [x; feature_tmp];
            end
        end
        
        % find pts index with valid value of each feature at this time point
        % for each pt, if one of the 71 features contains NaN, the pt will not be used as dataset for the model of this time range
        index = find(~isnan(mean(x)));
        pts_indicators = [pts_indicators,index];
        X = x(:,index)';
        Y = dci_labels(index);
        
        % concate pts with valid value of each time point of this time range 
        X_all = [X_all;X];
        Y_all = [Y_all;Y];
    end
    X_all_range{i,1} = X_all;
    Y_all_range{i,1} = Y_all;
    pts_indicators_range{i,1} = pts_indicators;
end

%%
num_pts_unique = cell(length(pts_indicators_range),1);
for i=1:33
    pts_indicators = pts_indicators_range{i};
    pts_unique = unique(pts_indicators);
    num_pts_unique{i} = length(pts_unique);
end
%%
% mypool = parpool(4);
% paroptions = statset('UseParallel',true);
%

%%

for rf_para = 50%40:10:60%50:10:100 %60
    for itrials = 1:num_itrials
        preds = cell(length(time_range),num_fold);
        labels = cell(length(time_range),num_fold);
        pred_probability = cell(length(time_range),num_fold);
        label_ratio = cell(length(time_range),num_fold);
        for i = 1:length(time_range)
            pts_indicators = pts_indicators_range{i};
            pts_unique = unique(pts_indicators);
            num_pts = length(pts_unique);
            num_test = round(num_pts/num_fold);
            idx = randperm(num_pts);
            pts_unique = pts_unique(idx);
        
            X_all = normalize(X_all_range{i});% normalize each column/feature
            Y_all = Y_all_range{i};
            
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
                test_label = Y_all_range{i}(test_index);
                train_data = X_all(train_index,:);
                train_label = Y_all_range{i}(train_index);
                ratio = sum(train_label)/length(train_label);
                train_data_ct = length(train_label);
                ct_ratio = train_data_ct/length(test_label);

                %% Random Forest
                B=TreeBagger(rf_para,train_data,train_label,'Method','Classification','OOBVarImp','On');
                weight_tmp=B.OOBPermutedVarDeltaError;
                [pred, probabilities] = predict(B, test_data);

                preds{i,ifold} = pred;
                labels{i,ifold} = test_label;
                pred_probability{i,ifold} = probabilities;
                label_ratio{i, ifold} = ratio;
                train_data_count{i, ifold} = train_data_ct;
                train_test_ratio{i, ifold} = ct_ratio;
            end
        end
        preds_all{itrials} = preds;
        labels_all{itrials} = labels;
        pred_probability_all{itrials} = pred_probability;

        dataset_info_all{itrials,1} = label_ratio;
        dataset_info_all{itrials,2} = train_data_count;
        dataset_info_all{itrials,3} = train_test_ratio;
    end
%%
    save(['E:\Zhongwei\SAH Code Publish\models\Baseline_RF\rf_',num2str(rf_para),'_trial',num2str(num_itrials), ...
        '_seed',num2str(seed),'_allFea'],'preds_all','labels_all','pred_probability_all','time_range');
end

