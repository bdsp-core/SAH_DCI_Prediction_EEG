
clear
close all
% fix features after DCI with features at DCI
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\vascular_features.mat');
pts_info = readtable('sah_dci_FINALadjud.csv');
pts_info.SID_str = arrayfun(@(x) sprintf('%03d', x), pts_info.SID, 'UniformOutput', false);
time_shift_all = nan(length(unique_names),1);
for ipatient = 1:length(unique_names)
    if ipatient==108
        continue
    end
    indx = contains(pts_info.SID_str,unique_names{ipatient});
    indx_tmp = find(indx==true);
    dates = datetime(char(pts_info.DayOfBleed(indx_tmp(1))),'Format','MM/dd/yyyy');
    times = datenum(char(pts_info.TimeOfBleed(indx_tmp(1))));
    dates = datetime(dates,'Format','MM/dd/yyyy HH:mm:SS');
    times = datetime(times,'ConvertFrom','datenum','Format','MM/dd/yyyy HH:mm:SS');
    timeBleed = dates+timeofday(times);
    
    if dci_labels(ipatient)==1
        dates = datetime(char(pts_info.Final_DCIDate(indx_tmp(1))),'Format','MM/dd/yyyy');
        times = datenum(char(pts_info.Final_DCITime(indx_tmp(1))));
        dates = datetime(dates,'Format','MM/dd/yyyy HH:mm:SS');
        times = datetime(times,'ConvertFrom','datenum','Format','MM/dd/yyyy HH:mm:SS');
        timeDCI = dates+timeofday(times);
        [h,m,s]= hms(time(between(timeBleed,timeDCI,'time')));
        time_shift_all(ipatient) = h;
    end
end


vascular_features = features;
% interpolate missing data
for idim = 1:size(vascular_features,1)
    for jdim = 1:size(vascular_features,2)
        feature_all = vascular_features{idim,jdim};
        for i = 1:size(feature_all,1)
            t = find(~isnan(feature_all(i,:)));
            yt = feature_all(i,t);
            f = fit(t',yt','linearinterp');
            y = feval(f,1:size(feature_all,2));
            y(y<prctile(yt,10)) = prctile(yt,10);
            y(y>prctile(yt,90)) = prctile(yt,90);
            if ~isnan(time_shift_all(i))
                y(time_shift_all(i):end) = mean(y(1:time_shift_all(i))); % mean?
            end
            feature_all(i,:) = y;
        end
        vascular_features{idim,jdim} = feature_all;
    end
end

load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\asymmetry_subtraction_features.mat')
% interpolate missing data
for idim = 1:size(asym_features,1)
    for jdim = 1:size(asym_features,2)
        feature_all = asym_features{idim,jdim};
        for i = 1:size(feature_all,1)
            t = find(~isnan(feature_all(i,:)));
            yt = feature_all(i,t);
            f = fit(t',yt','linearinterp');
            y = feval(f,1:size(feature_all,2));
            y(y<prctile(yt,10)) = prctile(yt,10);
            y(y>prctile(yt,90)) = prctile(yt,90);asym_features
            if ~isnan(time_shift_all(i))
                y(time_shift_all(i):end) = mean(y(1:time_shift_all(i)));
            end
            feature_all(i,:) = y;
        end
        asym_features{idim,jdim} = feature_all;
    end
end

load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\all_features_NEWspike30.mat')
% interpolate missing data
for idim = 1:size(features,1)
    for jdim = 1:size(features,2)
        feature_all = features{idim,jdim};
        for i = 1:size(feature_all,1)
            t = find(~isnan(feature_all(i,:)));
            yt = feature_all(i,t);
            f = fit(t',yt','linearinterp');
            y = feval(f,1:size(feature_all,2));
            y(y<prctile(yt,10)) = prctile(yt,10);
            y(y>prctile(yt,90)) = prctile(yt,90);
            if ~isnan(time_shift_all(i))
                y(time_shift_all(i):end) = mean(y(1:time_shift_all(i)));
            end
            feature_all(i,:) = y;
        end
        features{idim,jdim} = feature_all;
    end
end
feature_all = spike;
for i = 1:size(feature_all,1)
    t = find(~isnan(feature_all(i,:)));
    yt = feature_all(i,t);
    f = fit(t',yt','linearinterp');
    y = feval(f,1:size(feature_all,2));
    y(y<prctile(yt,10)) = prctile(yt,10);
    y(y>prctile(yt,90)) = prctile(yt,90);
    if ~isnan(time_shift_all(i))
        y(time_shift_all(i):end) = mean(y(1:time_shift_all(i)));
    end
    feature_all(i,:) = y;
end
spike = feature_all;

time_step = 1;
time_range = 43:time_step:240;

feature_all = cell(1, length(time_range));
label_all = cell(1, length(time_range));
auc_individual_feature = [];
for i = 1:length(time_range)
 
    iblock = 1;
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
                feature_tmp = asym_features{ifea,jfea}(:,time_range(i)-iblock+1)';
                x = [x; feature_tmp];
            end
        end
        for ifea = 4 % asymmetry delta mca
            for jfea = [2]
                feature_tmp = asym_features{ifea,jfea}(:,time_range(i)-iblock+1)';
                x = [x; feature_tmp];
            end
        end
        for ifea = 1 % asymmetry shannon aca mca
            for jfea = [1 2]
                feature_tmp = asym_features{ifea,jfea}(:,time_range(i)-iblock+1)';
                x = [x; feature_tmp];
            end
        end
        for ifea = [3 5] % asymmetry theta total mca
            for jfea = [2]
                feature_tmp = asym_features{ifea,jfea}(:,time_range(i)-iblock+1)';
                x = [x; feature_tmp];
            end
        end
        %% cumulative features adr pav spikes
        feature_tmp = cumsum(spike(:,time_range(1)-time_step:time_range(i)-iblock+1),2);
        x = [x; feature_tmp(:,end)'];
        feature_tmp = cumsum(features{2}(:,time_range(1)-time_step:time_range(i)-iblock+1),2);
        x = [x; feature_tmp(:,end)'];
        feature_tmp = cumsum(features{7}(:,time_range(1)-time_step:time_range(i)-iblock+1),2);
        x = [x; feature_tmp(:,end)'];
        
        
        index = find(~isnan(mean(x)));
        X = x(:,index)';
        Y = dci_labels(index);
        
        feature_all{i} = X;
        label_all{i} = Y;
end

%%
num_itrials = 1;
time_step = 6;
time_range = 48:time_step:240;
num_fold = 5;
seed = 19284;
rng(seed);

for rf_para = 50%20:10:100
for itrials = 1:num_itrials

preds = cell(length(time_range),num_fold);
labels = cell(length(time_range),num_fold);
pred_probability = cell(length(time_range),num_fold);
num_pts = 113;
num_test = round(num_pts/num_fold);
pts_unique = 1:num_pts;
idx = randperm(num_pts);
pts_unique = pts_unique(idx);
pts_fold = {};
for ifold = 1:num_fold
    if ifold~=num_fold
        start_index = (ifold-1)*num_test+1;
        end_index = ifold*num_test;
    else
        start_index = (ifold-1)*num_test+1;
        end_index = num_pts;
    end
    pts_fold{ifold} = pts_unique(start_index:end_index);
    
    test_index = pts_fold{ifold};
    train_index = setdiff(1:num_pts,test_index);
    
    auc_per_feature = zeros(size(feature_all{1},2), size(feature_all,2));
    diff_per_feature = zeros(size(feature_all{1},2), size(feature_all,2));
    best_time_all = [];
    for i = 1:length(time_range)
        X_all = [];
        Y_all = [];
        best_auc_all = [];
        best_diff_all = [];
        for iblock = time_step:-1:1
            X = feature_all{i*time_step-iblock+1};
            Y = label_all{i*time_step-iblock+1};
            

            
            % optimize each feature based on marginal seperation
            
            for ifea = 1:size(X,2)
                x_tmp = X(train_index,ifea);
                y_tmp = Y(train_index);
                % difference between two groups of a feature of each hour
                diff_per_feature(ifea,i*time_step-iblock+1) = (abs(mean(x_tmp(y_tmp==1))-mean(x_tmp(y_tmp==0))))/sqrt((std(x_tmp(y_tmp==1))*std(x_tmp(y_tmp==0))));
            end

            % max difference of hour_1 to hour_current
            [best_diff,best_time] = max(diff_per_feature');
            best_time_all = [best_time_all;best_time];
            % optimize each feature based on diff
            % keep best feature from hour_1 to hour_current
            x = [];
            for itime = 1:length(best_time)
                x = [x feature_all{best_time(itime)}(:,itime)];
                
            end
            X_all = [X_all;x];
            Y_all = [Y_all;label_all{best_time(itime)}];
            best_diff_all = [best_diff_all best_diff'];
            
        end
        X_all = normalize(X_all);
        test_index_all = [];
        train_index_all = [];
        for ipt = 1:time_step
            test_index_all = [test_index_all, test_index+num_pts*(ipt-1)];
            train_index_all = [train_index_all, train_index+num_pts*(ipt-1)];
        end
        
        test_data = X_all(test_index_all,:);
        test_label = Y_all(test_index_all);
        train_data = X_all(train_index_all,:);
        train_label = Y_all(train_index_all);

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
pts_fold_all{itrials} = pts_fold;
end
save(['E:\Zhongwei\SAH Code Publish\models\CarryForward\','cf_',num2str(rf_para),'_trial',num2str(num_itrials),'_seed',num2str(seed),'_slect_interp_cum'],'preds_all','labels_all','pred_probability_all','time_range','pts_fold');
end

%%
