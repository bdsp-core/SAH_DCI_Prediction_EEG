%Risk calibration curves of models at different time intervals after SAH

clear
load('E:\Zhongwei\SAH Code Publish\models\CarryForward\cf_50_trial_1_seed_19284_slect_interp_cum_.mat');
time_step = 6;
time_range = 48:time_step:240;
num_fold = 5;
%% plot calibration curve for hour48-72
pred_all = preds_all;

for itrials = 1:length(preds_all)
    labels = labels_all{itrials};
    pred_probability = pred_probability_all{itrials};
    preds = pred_all{itrials};
    label_all = [];
    pred_prob_all = [];
    preds_all =[];
    for i = 1:length(time_range)
        label_tmp = [];
        pred_prob_tmp = [];
        preds_tmp = [];
    
        for j = 1:num_fold
            num_pts = size(pred_probability{i,j},1)/time_step;
            for k = 1:num_pts
            label_tmp = [label_tmp;labels{i,j}(k:num_pts:end,:)'];
            pred_prob_tmp = [pred_prob_tmp;pred_probability{i,j}(k:num_pts:end,2)'];
            preds_tmp = [preds_tmp;preds{i,j}(k:num_pts:end,:)'];
            end
        end
        label_all = [label_all,label_tmp];
        pred_prob_all = [pred_prob_all,pred_prob_tmp];
        preds_all = [preds_all,preds_tmp];
        
    end
end
% plot calibration curve at 48-72 h after DCI
num = 1;
Brier= [];
M = 3;
% M = 8;
xx = linspace(0, 1, M+1); 
figure,
plot(xx, xx, 'k--')
hold on
markerlist = {'r-o','g-+','b-^','m-s','c-d'};
for plot_time = 48:6:72 
    
    y = label_all(:,plot_time-43-4:plot_time-43+1);
    yp = pred_prob_all(:,plot_time-43-4:plot_time-43+1);
    
    y = reshape(y,size(y,1)*size(y,2),1);
    yp = reshape(yp,size(yp,1)*size(yp,2),1);
    
    % Compute Calibarition: Fraction of positive vs. prediction %
    
    xx = linspace(0, 1, M+1); 
    yy = NaN(M,1);     
    for k = 1:M
        y_k = y(yp>=xx(k) & yp<xx(k+1));
        yy(k) = sum(y_k)/length(y_k);
    end

    xx = linspace(0, 1, M+1);
    xx = (xx(1:M)+xx(2:end))/2;
    Brier(num) = sum(abs(xx-yy'))/length(xx); %% Brier's Score
    eval(['h',num2str(num),'=plot(xx, yy,markerlist{num},''linewidth'',2);']);
    axis square
    box on
    grid on
    title(['Risk Calibration'])
    xlabel('Prediction')
    ylabel('Proportion of the Positive')
    
    num = num+1;
end
legend([h1 h2 h3 h4 h5],{['48h,',num2str(Brier(1))],['54h',num2str(Brier(2))],['60h',num2str(Brier(3))],['66h',num2str(Brier(4))],['72h',num2str(Brier(5))]});

%% plot calibration curve 12 h prior to DCI time

pts_index = cat(2,pts_fold{:})';
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

y = [];
yp = [];
for i = 1:length(time_shift_all)
    if ~isnan(time_shift_all(i))
        if (time_shift_all(i)-12)<=size(pred_prob_all,2)&&time_shift_all(i)-12-5>=1
            yp = [yp; pred_prob_all(i,time_shift_all(i)-12-5:time_shift_all(i)-12)']; %use probability 12 h prior to dci time
            y = [y; label_all(i,time_shift_all(i)-12-5:time_shift_all(i)-12)'];
        else if (time_shift_all(i)-12)>size(pred_prob_all,2)
            yp = [yp; pred_prob_all(i,end-5:end)'];
            y = [y; label_all(i,end-5:end)'];
        else
            yp = [yp; pred_prob_all(i,1:time_shift_all(i)-12)'];
            y = [y; label_all(i,1:time_shift_all(i)-12)'];
            end
        end
    else
        yp = [yp; pred_prob_all(i,end-5:end)']; % use probability at 240h for no dci pts
        y = [y; label_all(i,end-5:end)'];
    end
end

% Compute Calibration: Fraction of positive vs. Predition %
M = 3;
xx = linspace(0, 1, M+1); 
yy = NaN(M,1);     
for k = 1:M
    y_k = y(yp>=xx(k) & yp<xx(k+1));
    yy(k) = sum(y_k)/length(y_k);
end

figure;
hold on
    plot(xx, xx, 'k--')
    xx = linspace(0, 1, M+1);
    xx = (xx(1:M)+xx(2:end))/2;
Brier = sum(abs(xx-yy'))/length(xx); %% Brier's Score
    h1 = plot(xx, yy,'b-^','linewidth',2);
    
    axis square
    box on
    grid on
    title('Risk Calibration at 18h-12h Prior to DCI')
    xlabel('Prediction')
    ylabel('Proportion of the Positive')
    legend(h1,['18h-12h Prior to DCI,',num2str(Brier)]);
hold off