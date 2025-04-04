% Fig 2 abc
%% load prepared features
load('./ExtractFeatures/all_features_NEWspike30.mat','dci_labels','unique_names_dci','features','spike');

ifeature = 2;
feature_all = features{ifeature};
num_pts = size(feature_all,1);

feature_names = {'Shannon'
    'AlphaDeltaRatio'
    'TotalPow'
    'DeltaPow'
    'ThetaPow'
    'AlphaPow'
    'PerAlphaVar'
    };

%%
% feature_all = 10*log10(feature_all);
index_tmp = find(dci_labels==0);
feature_good = feature_all(dci_labels==0,:);
feature_mean = nanmean(feature_good');
[a,index_good]=sort(feature_mean);
index_good = index_tmp(index_good);
feature_good = feature_all(index_good,:);
unique_names_good = unique_names_dci(index_good);

index_tmp = find(dci_labels==1);
feature_bad = feature_all(dci_labels==1,:);
feature_mean = nanmean(feature_bad');
[a,index_bad]=sort(feature_mean);
index_bad = index_tmp(index_bad);
feature_bad = feature_all(index_bad,:);
unique_names_bad = unique_names_dci(index_bad);

unique_names1 = [unique_names_good';unique_names_bad'];
feature_all1 = [feature_good;feature_bad];

range_start = 24;
range_end = 10*24+1;


%% plot distribution of two groups
% Plot No-DCI Feature Trend
for i = 1:size(feature_good,2)
    feature_value = feature_good(:,i);
    feature_value(feature_value>prctile(feature_value,90)) = prctile(feature_value,90);
    mean_feature(i) = nanmean(feature_value);
    median_feature(i) = nanmedian(feature_value);
    std_feature(i) = nanstd(feature_value);
    valid_pts(i) = sum(~isnan(feature_value));
end

x= range_start:range_end;
y = mean_feature(range_start:range_end);
xb = x; yb = y; 
f = fit(x',y','smoothingspline','SmoothingParam',0.001);
figure,h1 = plot(f,'b');
hold on,
scatter(x,y,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)


% Compute Mean and STD for DCI Group
for i = 1:size(feature_bad,2)
    feature_value = feature_bad(:,i);
    feature_value(feature_value>prctile(feature_value,90)) = prctile(feature_value,90);
    mean_feature(i) = nanmean(feature_value);
    std_feature(i) = nanstd(feature_value);
    valid_pts(i) = sum(~isnan(feature_value));
end

% plot p values indicators 
% statistical Significance Testing (t-test): Independent t-test (ttest2) compares No-DCI vs. DCI feature values for each time point.
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

% Filter Out Isolated Significant Points
% Count p-value as “significant” only if p<0.05
% A point is considered significant only if at least 2 consecutive time points show significance.
indx = find(p_indicator==1);
a = find(diff(indx)==1);
valid_indx = unique([indx(a),indx(a)+1]);
p_indicator_new = p_indicator;
p_indicator_new(setdiff(indx,valid_indx)) = 0;

x= range_start:range_end;
y = mean_feature(range_start:range_end);
xr=x; yr = y; 
f = fit(x',y','smoothingspline','SmoothingParam',0.001);
hold on,h2 = plot(f,'r');hold on,
scatter(x,y,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

% Mark Significant Time Points
p_indicator(p_indicator==0) = nan;
p_indicator_new(p_indicator_new==0) = nan;
h3 = plot(range_start:range_end,p_indicator_new(range_start:range_end)*0.02,'k-','LineWidth',2);

% scatter(range_start:range_end,p_indicator(range_start:range_end)*300,'k*');
xlim([range_start-1,range_end]),
xlabel('Days After Bleed')
set(gca, 'Xtick',range_start:24:range_end,'Xticklabel',num2cell(1:10))
ylabel(feature_names{ifeature});
legend([h1 h2 h3],{'No DCI','DCI','p<0.05'});
legend boxoff


%% plot with smoother curve fit(fig in paper)
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
h3 = plot(xr,h*0.025,'k','linewidth',2)

set(gca, 'Xtick',range_start:24:range_end,'Xticklabel',num2cell(1:10))
% ylabel('Spike Frequency')
ylabel(feature_names{ifeature});
xlabel('Days After Bleed')

legend([h1 h2 h3],{'No DCI','DCI','p < 0.05'},'Location','northeast');
legend boxoff
set(gca,'tickdir','out')

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
