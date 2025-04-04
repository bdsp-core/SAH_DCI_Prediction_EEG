% Plot heatmap and line plot for cumulative features (PAV in this example)
% Fig - 2.d

%% load prepared features
clear
% close all;
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\all_features_NEWspike30.mat','dci_labels','unique_names_dci','features');
ifeature = 7;
feature_all = features{ifeature};

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


% trunc data of 10 days
range_start = 24;
range_end = 10*24;
feature_all = cumsum(feature_all(:,range_start:range_end),2);

dci_labels = dci_labels;
unique_names_dci = unique_names_dci;
num_pts = size(feature_all,1);

feature_names = {'Shannon'
    'AlphaDeltaRatio'
    'TotalPow'
    'DeltaPow'
    'ThetaPow'
    'AlphaPow'
    'PerAlphaVar'
    };


% feature good/bad set up for showing two class on half of the pic
% feature_all = 10*log10(feature_all);
index_tmp = find(dci_labels==0);
feature_good = feature_all(dci_labels==0,:);
feature_mean = nanmean(feature_good');
[a,index_good]=sort(feature_mean);
index_good = index_tmp(index_good);
feature_good = feature_all(index_good,:);
unique_names_good = unique_names_dci(index_good);
% score_good = score_all(index_good,:);

index_tmp = find(dci_labels==1);
feature_bad = feature_all(dci_labels==1,:);
feature_mean = nanmean(feature_bad');
[a,index_bad]=sort(feature_mean);
index_bad = index_tmp(index_bad);
feature_bad = feature_all(index_bad,:);
unique_names_bad = unique_names_dci(index_bad);
% score_bad = score_all(index_bad,:);

% feature_bad = feature_bad(1:end-3,:);
unique_names1 = [unique_names_good';unique_names_bad'];
feature_all1 = [feature_good;feature_bad];
% score_all = [score_good;score_bad];

range_start = 1;
range_end = size(feature_all1,2);

%% heatmap representation

fig = figure('units','normalized','outerposition',[0 0 1 1]);
mn = prctile(feature_all1(:),10); mx = prctile(feature_all1(:),90);% 90
imagesc(feature_all1,[mn,mx]);
% colormap(flipud(cold));
set(gca, 'YTick', 1:num_pts, 'YTickLabel', unique_names1);
title('PAV Frequency')
% set(gca, 'YTick', (1/num_pts:1/num_pts:1)*num_pts, 'YTickLabel', unique_names1)
xlabel('Hours After Bleed');
ylabel('Patients/Subjects');
% title(feature_names{ifeature});
colorbar;

xlim([range_start,range_end]),xlabel('Days After Bleed')
set(gca, 'Xtick',range_start:24:range_end,'Xticklabel',num2cell(1:10))


%% Line Plots for Mean and Variability Analysis
% visualizes mean feature values with standard error bars for DCI and No-DCI groups
figure

% mean feature values with standard error bars for DCI
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


% mean feature values with standard error bars for non-DCI
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

xlim([range_start-1,range_end]),xlabel('Days After Bleed')
set(gca, 'Xtick',range_start:24:range_end,'Xticklabel',num2cell(1:10))
ylim([0 30])
ylabel('Feature Values')
title('PAV Frequency');
axis tight;
legend([h1 h2],{'No DCI','DCI'});