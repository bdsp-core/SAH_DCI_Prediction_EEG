%% AUC for DCI “alarm” of identifying PRIOR to DCI (False+,False-,True+,True-)
clear
results = load('E:\Zhongwei\SAH Code Publish\models\CarryForward\cf_50_trial_1_seed_19284_slect_interp_cum_.mat');
time_step = 6;
time_range = 48:time_step:240;
num_fold = 5;

X = [];
Y = [];
AUC = [];
mean_curve = [];
mean_curve_seperate = [];
mean_all_trials = [];
mean_trials = [];
boostrap_times = 5000;
%% rearrange output matrix and get pt DCI onset time
% re-structure pred_prob_all and preds_all to (113,198) in order to later
% match the DCI onset time using pts_fold

% get time_shift for patients matched with the output
pts_index = cat(2,results.pts_fold{:})';
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\all_features_NEWspike30.mat','unique_names','dci_labels','unique_names_dci');
pts_info = readtable('.\sah_dci_FINALadjud.csv');
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
time_shift_all = time_shift_all(pts_index);
time_shift_all = time_shift_all-43; % start at hour 43 after SAH
dci_labels = dci_labels(pts_index);
unique_names_dci = unique_names_dci(pts_index);


%% Calculate AUC
for itrials = 1:length(results.preds_all)
    labels = results.labels_all{itrials};
    pred_probability = results.pred_probability_all{itrials};
    preds = results.preds_all{itrials};
    pred_prob_all = [];
    preds_all =[];
    pred_prob_ave = [];
    pred_prob_cum = [];

    for i = 1:length(time_range)
        pred_prob_tmp = [];
        preds_tmp = [];
    
        for j = 1:num_fold
            num_pts = size(pred_probability{i,j},1)/time_step;
            for k = 1:num_pts
            pred_prob_tmp = [pred_prob_tmp;pred_probability{i,j}(k:num_pts:end,2)'];
            preds_tmp = [preds_tmp;preds{i,j}(k:num_pts:end,:)'];
            end
        end

        pred_prob_all= [pred_prob_all,pred_prob_tmp];
        preds_all= [preds_all,preds_tmp];
        
        % average past predicted proba
        if i==1
            pred_prob_cum = pred_prob_tmp;
        else
            pred_prob_cum = pred_prob_cum + pred_prob_tmp;
        end

        pred_prob_ave = [pred_prob_ave, pred_prob_cum./i];

    end

    % bootstrap
    for iboost = 1:boostrap_times
        [pred_prob_all_boostrap,idx] = datasample(pred_prob_ave,size(pred_prob_ave,1),1);
        unique_names_boostrap = unique_names(idx);
        time_shift_all_boostrap = time_shift_all(idx);
        dci_labels_boostrap = dci_labels(idx);
    
        TPR = [];
        FPR = [];
        for threshold = 0.05:0.05:0.95
            TP_count = 0;
            FP_count = 0;
            for i = 1:length(unique_names_boostrap)
                prediction = pred_prob_all_boostrap(i,:)>threshold;
                % only count prediction of DCI patient before DCI onset
                if ~isnan(time_shift_all_boostrap(i))
                    if time_shift_all_boostrap(i)<=198
                        prediction = prediction(1:time_shift_all_boostrap(i));
                    else
                        prediction = prediction(1:end);
                    end
                    if ismember(1,prediction)
                        TP_count = TP_count+1;
                    end
                else
                    if ismember(1,prediction)
                        FP_count = FP_count+1;
                    end
                end
                    
            end
            TPR = [TPR;TP_count/(sum(dci_labels_boostrap==1))];
            FPR = [FPR;FP_count/(sum(dci_labels_boostrap==0))];
    
        end
        [a,b] = sort(FPR);
        FPR = FPR(b);
        TPR = TPR(b);
    
        X{iboost} = [0;TPR;1];
        Y{iboost} = [0;FPR;1];
        AUC{iboost} = trapz([0;FPR;1],[0;TPR;1]);
    
        x_adj = adjust_unique_points([0;FPR;1]); %interp1 requires unique points
        if iboost ==1
            mean_curve = interp1(x_adj, [0;TPR;1], linspace(0, 1, 100));
        else
            mean_curve = mean_curve+ (interp1(x_adj, [0;TPR;1], linspace(0, 1, 100)));
        end
        mean_curve_seperate = [mean_curve_seperate;(interp1(x_adj, [0;TPR;1], linspace(0, 1, 100)))];
    end

    % mean result of all trials
    if itrials==1
        mean_all_trials = mean_curve;
    else
        mean_all_trials = mean_all_trials + mean_curve;
    end
    % seperately store results of different trials
    mean_trials = [mean_trials; mean_curve];
end

%%
yCI95_low = prctile(mean_curve_seperate,2.5);
yCI95_high = prctile(mean_curve_seperate,97.5);
yCI95 = [mean_curve/boostrap_times-yCI95_low;yCI95_high-mean_curve/boostrap_times]';

%%
% ROC curve max carry forward model considering only data up through the time of DCI to represent prediction performance prior to DCI onset. Shaded areas denote 95 % confidence intervals.
figure,
[l4,p] = boundedline(linspace(0, 1, 100), mean_curve/boostrap_times, yCI95,'-r','alpha');% --b* --ro --g+ --cs

hold on
plot(linspace(0, 1, 100), linspace(0, 1, 100), 'k--')
axis square
box on
title('Prediction Performance Considering DCI time - average past proba')
xlabel('1-Specificity')
ylabel('Sensitivity')
text(0.6,0.2,['AUC = ',num2str(mean(cell2mat(AUC)))],'FontSize',12);
text(0.6,0.1,['(95% CI:',num2str(prctile(cell2mat(AUC),2.5)),'-',num2str(prctile(cell2mat(AUC),97.5))],'FontSize',12);
hold off
            
%%
TPR = [];
FPR = [];
for threshold = 0.05:0.05:0.95
TP_count = 0;
FP_count = 0;
for i = 1:length(unique_names)
    prediction = pred_prob_ave(i,:)>threshold;
    if ~isnan(time_shift_all(i))
        if time_shift_all(i)<=198
            prediction = prediction(1:time_shift_all(i));
        else
            prediction = prediction(1:end);
        end
        if ismember(1,prediction)
            TP_count = TP_count+1;
        end
    else
        if ismember(1,prediction)
            FP_count = FP_count+1;
        end
    end
        
end
TPR = [TPR;TP_count/(sum(dci_labels==1))];
FPR = [FPR;FP_count/(sum(dci_labels==0))];

end
[a,b] = sort(FPR);
FPR = FPR(b);
TPR = TPR(b);

X = [0;TPR;1];
Y = [0;FPR;1];
AUC = trapz(Y,X);
figure
plot(Y,X,'b-','linewidth',2)
hold on,plot([0;1],[0;1],'--k')
axis square
    box on
    grid on
xlabel('False positive rate');ylabel('True positive rate');
title('Prediction Performance Considering DCI time - average past proba')
text(0.6,0.1,['AUC = ',num2str(AUC)],'FontSize',12);


%%
function x= adjust_unique_points(Xroc)
    x= zeros(1, length(Xroc));
    aux= 0.0001;
    for i=1: length(Xroc)
        if i~=1
            x(i)= Xroc(i)+aux;
            aux= aux+0.0001;
        end
        
    end
end