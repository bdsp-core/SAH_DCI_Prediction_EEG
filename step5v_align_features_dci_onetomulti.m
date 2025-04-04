%% Visualization
% heatmap and ttest of aligned features and cumulative features.

clear
close all;
% change title, ifeature index
% for vascular and asymmetry:
% load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_asym.mat');
% % features = vas_features_dci;
% features = asym_features_dci;
% ifeature = 2;
% ivascular = 1;
% feature_all = features{ifeature,ivascular};

% for global features:
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_global.mat');
ifeature = 2;
feature_all = features_dci{ifeature};

load('ctrl_matchtime_dci.mat');
total_cases_after_dci = length(cat(1,match_ctrls{:}))+length(ctrl_sids);
dci_labels = dci_labels_new;

pts_info = readtable('.\sah_dci_FINALadjud.csv');
pts_info.SID_str = arrayfun(@(x) sprintf('%03d', x), pts_info.SID, 'UniformOutput', false);
time_shift_all = [];
dci_labels_new = [];
pts_index = [];
ctrl_index = find(dci_labels==0);
unique_names_dci_all = {};
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
        time_shift_all = [time_shift_all,h];
        dci_labels_new = [dci_labels_new,1];
        pts_index = [pts_index, ipatient];
        unique_names_dci_all = [unique_names_dci_all,unique_names_dci{ipatient}];
    end
end

for i = 1:length(match_ctrls)
    for j = 1:length(match_ctrls{i})

        pts_index = [pts_index, ctrl_index(match_ctrls{i}(j))];
        time_shift_all = [time_shift_all,time_shift_all(i)];
        dci_labels_new = [dci_labels_new,0];
        unique_names_dci_all = [unique_names_dci_all,unique_names_dci{ctrl_index(match_ctrls{i}(j))}];

    end
end

dci_point = 1000; % each pt will be aligned to dci onset time point=1000
num_pts = size(feature_all,1);

feature_names = {'Shannon'
    'AlphaDeltaRatio'
    'TotalPow'
    'DeltaPow'
    'ThetaPow'
    'AlphaPow'
    'PerAlphaVar'
    };
vascular_names = {'ACA\_left'
    'MCA\_left'
    'PCA\_left'
    'ACA\_right'
    'MCA\_right'
    'PCA\_right'
    };


% feature_all = 10*log10(feature_all);
index_tmp = find(dci_labels_new==0);
feature_good = feature_all(dci_labels_new==0,:);
feature_mean = nanmean(feature_good');
[a,index_good] = sort(feature_mean);
index_good = index_tmp(index_good);
feature_good = feature_all(index_good,:);
unique_names_good = unique_names_dci_all(index_good);
% score_good = score_all(index_good,:);

index_tmp = find(dci_labels_new==1);
feature_bad = feature_all(dci_labels_new==1,:);
feature_mean = nanmean(feature_bad');
[a,index_bad] = sort(feature_mean);
index_bad = index_tmp(index_bad);
feature_bad = feature_all(index_bad,:);
unique_names_bad = unique_names_dci_all(index_bad);
% score_bad = score_all(index_bad,:);

unique_names1 = [unique_names_good';unique_names_bad'];
feature_all1 = [feature_good;feature_bad];
% score_all = [score_good;score_bad];

range_start = 1000-5*24;range_end = 1000+3*24;

fig = figure('units','normalized','outerposition',[0 0 1 1]);
mn = prctile(feature_all1(:),10); mx = prctile(feature_all1(:),90);% 90
imagesc(feature_all1(:,range_start:range_end),[mn,mx]);
colormap;
set(gca, 'YTick', (1/num_pts:1/num_pts:1)*num_pts, 'YTickLabel', unique_names1)
set(gca, 'XTick', 1:24:size(feature_all1,2), 'XTickLabel', cellstr(num2str((((range_start:24:range_end)-dci_point)/24)')));
xlabel('Days After DCI');
% title([feature_names{ifeature},'\_',vascular_names{ivascular}]);
title(feature_names{ifeature});
% title('Spike Frequency')
colorbar;

%% plot distribution of two groups
% figure
for i = 1:size(feature_good,2)
    feature_value = feature_good(:,i);
    feature_value = feature_value(feature_value>=prctile(feature_value,10)&feature_value<=prctile(feature_value,90));
    mean_feature(i) = nanmean(feature_value);
    std_feature(i) = nanstd(feature_value);
    valid_pts(i) = sum(~isnan(feature_value));
end

x= range_start:range_end;
y = mean_feature(range_start:range_end);
xb = x; yb = y; 
f = fit(xb',yb','smoothingspline','SmoothingParam',0.0001);
fig = figure('units','normalized','outerposition',[0 0 1 1]);
h1 = plot(f,'b');
hold on,
scatter(x,y,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

for i = 1:size(feature_bad,2)
    feature_value = feature_bad(:,i);
    feature_value = feature_value(feature_value>=prctile(feature_value,10)&feature_value<=prctile(feature_value,90));
    mean_feature(i) = nanmean(feature_value);
    std_feature(i) = nanstd(feature_value);
    valid_pts(i) = sum(~isnan(feature_value));
end

x= range_start:range_end;
y = mean_feature(range_start:range_end);
xr=x; yr = y; 
f = fit(x',y','smoothingspline','SmoothingParam',0.0001);
hold on,h2 = plot(f,'r');hold on,
scatter(x,y,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

% plot p values indicators
p_values = nan(1,size(feature_good,2));
p_indicator = nan(1,size(feature_good,2));
for itime = 1:size(feature_good,2)
    if sum(~isnan(feature_good(:,itime)))==0||sum(~isnan(feature_bad(:,itime)))==0
        continue
    else
        [h,p] = ttest2(feature_good(:,itime),feature_bad(:,itime));
        p_indicator(itime) = h;
        p_values(itime) = p;
    end
end
% Count p-value as “significant” only if p<0.05 and this is true for at least 2 consecutive hours
indx = find(p_indicator==1);
a = find(diff(indx)==1);
valid_indx = unique([indx(a),indx(a)+1]);
p_indicator_new = p_indicator;
p_indicator_new(setdiff(indx,valid_indx)) = 0;
% Mark Significant Time Points
p_indicator(p_indicator==0) = nan;
p_indicator_new(p_indicator_new==0) = nan;

plot(range_start:range_end,p_indicator_new(range_start:range_end)*0.035,'k-','LineWidth',8);

xlim([range_start-1,range_end]),xlabel('Days After DCI')
set(gca, 'XTick', range_start:24:range_end, 'XTickLabel', cellstr(num2str((((range_start:24:range_end)-dci_point)/24)')));
% ylabel(feature_names{ifeature})
title(feature_names{ifeature});
% title([feature_names{ifeature},'\_',vascular_names{ivascular}]);
% title('Spike Frequency');
% axis tight;
legend([h1 h2],{'No DCI','DCI'});
legend boxoff
% saveas(gcf,['D:\Research\SAH_DCI_Jenn\spectral_paper\AlignedToDCI\','Spike_Curve','.png']);
% saveas(gcf,['D:\Research\SAH_DCI_Jenn\results\alignWithDCI\', 'Spike','_Curve','.png']);
% close all;
%%
figure(2); clf; 
scatter(xb,yb,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
hold on;
scatter(xr,yr,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)


fb = fit(xb',yb','smoothingspline','SmoothingParam',0.0001);
h1=plot(fb,'b');
set(h1,'linewidth',2)

fr = fit(xr',yr','smoothingspline','SmoothingParam',0.0001);
h2=plot(fr,'r')
set(h2,'linewidth',2)
% ylim([0 300])
legend off
% p-values
w = 4; 
for i = 1:length(xr); 
    ind = find(abs(xr-xr(i))<=w); 
    xx = yr(ind); 
    yy = yb(ind); 
    [h(i),p(i)] = ttest(xx,yy)
end

ind = find(h==0); h(ind) = nan; 
h3 = plot(xr,h*25,'k','linewidth',2);

% set(gca, 'XTick', range_start:24:range_end, 'XTickLabel', cellstr(num2str((((range_start:24:range_end)-dci_point)/24)')));
set(gca, 'Xtick',range_start:24:range_end,'Xticklabel',num2cell(1:10))
% ylabel('Spike Frequency')
ylabel(feature_names{ifeature});
% xlim([range_start-1,range_end]),
xlabel('Days After DCI')
title(feature_names{ifeature});
% title([feature_names{ifeature},'\_',vascular_names{ivascular}]);
% title('Spike Frequency');

legend([h1 h2 h3],{'No DCI','DCI','p < 0.05'}); %,'Location','southeast'
legend boxoff
set(gca,'tickdir','out')








%% Visualization of cumulative features

% change title, ifeature index
% for vascular and asymmetry:
% load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_asym.mat');
% % features = vas_features_dci;
% features = asym_features_dci;
% ifeature = 2;
% ivascular = 1;
% feature_all = features{ifeature,ivascular};

% for global features:
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_global.mat');

ifeature = 2;%:7
% feature_all = features_dci{ifeature};
feature_all = spike_dci;
num_pts = size(feature_all,1);

%% interpolate missing data
for i = 1:size(feature_all,1)
    t = find(~isnan(feature_all(i,:)));
    yt = feature_all(i,t);
    f = fit(t',yt','linearinterp');
    y = feval(f,1:size(feature_all,2));
    y(y<prctile(yt,10)) = prctile(yt,10);
    y(y>prctile(yt,90)) = prctile(yt,90);
    feature_all(i,:) = y;
end

range_start = 1000-5*24;range_end = 1000+3*24;
% feature_all = feature_all(:,range_start:range_end);
feature_all = cumsum(feature_all(:,range_start:range_end),2);

feature_names = {'Shannon'
    'AlphaDeltaRatio'
    'TotalPow'
    'DeltaPow'
    'ThetaPow'
    'AlphaPow'
    'PerAlphaVar'
    };


% feature_all = 10*log10(feature_all);
index_tmp = find(dci_labels_new==0);
feature_good = feature_all(dci_labels_new==0,:);
feature_mean = nanmean(feature_good');
[a,index_good] = sort(feature_mean);
index_good = index_tmp(index_good);
feature_good = feature_all(index_good,:);
unique_names_good = unique_names_dci_all(index_good);
% score_good = score_all(index_good,:);

index_tmp = find(dci_labels_new==1);
feature_bad = feature_all(dci_labels_new==1,:);
feature_mean = nanmean(feature_bad');
[a,index_bad] = sort(feature_mean);
index_bad = index_tmp(index_bad);
feature_bad = feature_all(index_bad,:);
unique_names_bad = unique_names_dci_all(index_bad);
% score_bad = score_all(index_bad,:);

unique_names1 = [unique_names_good';unique_names_bad'];
feature_all1 = [feature_good;feature_bad];
% score_all = [score_good;score_bad];


range_start = 1;
range_end = size(feature_all1,2);

fig = figure('units','normalized','outerposition',[0 0 1 1]);
mn = prctile(feature_all1(:),10); mx = prctile(feature_all1(:),90);% 90
imagesc(feature_all1(:,range_start:range_end),[mn,mx]);
colormap();
set(gca, 'YTick', (1/num_pts:1/num_pts:1)*num_pts, 'YTickLabel', unique_names1)
set(gca, 'XTick', range_start:24:range_end, 'XTickLabel', cellstr(num2str([-5:1:3]')));
xlabel('Days After DCI');
title(feature_names{ifeature});
% title('Spike Frequency')
colorbar;

%%
figure
for i = 1:size(feature_good,2)
    feature_value = feature_good(:,i);
    mean_feature(i) = nanmean(feature_value);
    std_feature(i) = nanstd(feature_value);
    valid_pts(i) = sum(~isnan(feature_value));
end


[l,p] = boundedline(range_start:range_end, mean_feature(range_start:range_end), std_feature(range_start:range_end)./sqrt(valid_pts(range_start:range_end)),'--b*','alpha');% --b* --ro
h1 = outlinebounds(l,p);

mean_feature_good = mean_feature(range_start:range_end);
std_feature_good = std_feature(range_start:range_end);

for i = 1:size(feature_bad,2)
    feature_value = feature_bad(:,i);
    mean_feature(i) = nanmean(feature_value);
    std_feature(i) = nanstd(feature_value);
    valid_pts(i) = sum(~isnan(feature_value));
end

[l,p] = boundedline(range_start:range_end, mean_feature(range_start:range_end), std_feature(range_start:range_end)./sqrt(valid_pts(range_start:range_end)),'--ro','alpha');% --b* --ro
h2 = outlinebounds(l,p);

mean_feature_bad = mean_feature(range_start:range_end);
std_feature_bad = std_feature(range_start:range_end);

xlim([range_start-1,range_end]),xlabel('Days After DCI')
set(gca, 'XTick', range_start:24:range_end, 'XTickLabel', cellstr(num2str([-5:1:3]')));
ylabel('Feature Values')
title(feature_names{ifeature});
% title('Spike Frequency');
axis tight;
legend([h1 h2],{'No DCI','DCI'});
