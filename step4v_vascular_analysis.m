
clear
load('.\ExtractFeatures\vascular_features.mat')

%% plot heatmap of vascular territory of a feature
% of all patients in selected time range
ifeature = 6%:7;
ivascular = 6%1:6;
feature_all = features{ifeature,ivascular};
score_all = feature_score{ifeature,ivascular};
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
vascular_save = {'ACA_left'
    'MCA_left'
    'PCA_left'
    'ACA_right'
    'MCA_right'
    'PCA_right'
    };

feature_all = 10*log10(feature_all);
index_tmp = find(dci_labels==0);
feature_good = feature_all(dci_labels==0,:);
feature_mean = nanmean(feature_good');
[a,index_good]=sort(feature_mean);
index_good = index_tmp(index_good);
feature_good = feature_all(index_good,:);
unique_names_good = unique_names_dci(index_good);
score_good = score_all(index_good,:);

index_tmp = find(dci_labels==1);
feature_bad = feature_all(dci_labels==1,:);
feature_mean = nanmean(feature_bad');
[a,index_bad]=sort(feature_mean);
index_bad = index_tmp(index_bad);
feature_bad = feature_all(index_bad,:);
unique_names_bad = unique_names_dci(index_bad);
score_bad = score_all(index_bad,:);

unique_names1 = [unique_names_good';unique_names_bad'];
feature_all1 = [feature_good;feature_bad];
score_all = [score_good;score_bad];


% heatmap representaion 
fig = figure('units','normalized','outerposition',[0 0 1 1]);
mn = prctile(feature_all1(:),0); mx = prctile(feature_all1(:),100);% 90
imagesc(feature_all1,[mn,mx]);
colormap(flipud(parula));
colorbar;
set(gca, 'YTick', 1:num_pts, 'YTickLabel', unique_names1)
% set(gca, 'YTick', (1/num_pts:1/num_pts:1)*num_pts, 'YTickLabel', unique_names1)

title([feature_names{ifeature},'\_',vascular_names{ivascular}]);
xlabel('Hours After Bleed');
% title('Spike Frequency')
% caxis([-50,350])
xlim([40,200])

%% Plot Vascular territory feature discrimination
% SmoothingPara=0.001

% figure
mean_feature = nanmean(feature_good);
std_feature = nanstd(feature_good);
valid_pts = sum(~isnan(feature_good));
range_start = 24;
range_end = 10*24+1;
x= range_start:range_end;
y = mean_feature(range_start:range_end);
xb = x;
yb = y;
f = fit(x',y','smoothingspline','SmoothingParam',0.001);
fig = figure('units','normalized','outerposition',[0 0 1 1]);
h1 = plot(f,'b');
hold on,
scatter(x,y,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% [l,p] = boundedline(range_start:range_end, mean_feature(range_start:range_end), std_feature(range_start:range_end)./sqrt(valid_pts(range_start:range_end)),'--b*','alpha');% --b* --ro
% h1 = outlinebounds(l,p);

mean_feature = nanmean(feature_bad);
std_feature = nanstd(feature_bad);
valid_pts = sum(~isnan(feature_bad));
x= range_start:range_end;
y = mean_feature(range_start:range_end);
xr = x;
yr = y;
f = fit(x',y','smoothingspline','SmoothingParam',0.001);
hold on,h2 = plot(f,'r');hold on,
scatter(x,y,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% [l,p] = boundedline(range_start:range_end, mean_feature(range_start:range_end), std_feature(range_start:range_end)./sqrt(valid_pts(range_start:range_end)),'--ro','alpha');% --b* --ro
% h2 = outlinebounds(l,p);

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
p_indicator(p_indicator==0) = nan;
p_indicator_new(p_indicator_new==0) = nan;
plot(range_start:range_end,p_indicator_new(range_start:range_end)*8,'k-','LineWidth',8);

xlim([range_start,range_end]),xlabel('Days After Bleed')
set(gca, 'Xtick',range_start:24:range_end,'Xticklabel',num2cell(1:10))
ylabel([feature_names{ifeature},'\_',vascular_names{ivascular}])

legend([h1 h2],{'No DCI','DCI'});
legend boxoff
legend('Location','southeast')

%% plot smaller smoothing fit SmoothingPara=0.0001
fig = figure('units','normalized','outerposition',[0 0 1 1]);
% figure(2);
% clf; 
scatter(xb,yb,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
hold on;
scatter(xr,yr,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)


fb = fit(xb',yb','smoothingspline','SmoothingParam',0.0001);
h1=plot(fb,'b');
set(h1,'linewidth',2)

fr = fit(xr',yr','smoothingspline','SmoothingParam',0.0001);
h2=plot(fr,'r');
set(h2,'linewidth',2)
% ylim([0 300])
legend off
% p-values
w = 4; 
for i = 1:length(xr); 
    ind = find(abs(xr-xr(i))<=w); 
    xx = yr(ind); 
    yy = yb(ind); 
    [h(i),p(i)] = ttest(xx,yy);
end

ind = find(h==0); h(ind) = nan; 
h3 = plot(xr,h*1,'k','linewidth',2);

set(gca, 'Xtick',range_start:24:range_end,'Xticklabel',num2cell(1:10))
ylabel([feature_names{ifeature},'\_',vascular_names{ivascular}])
xlim([range_start,range_end]),xlabel('Days After Bleed')

legend([h1 h2 h3],{'No DCI','DCI','p < 0.05'},'Location','southeast');%
legend boxoff
set(gca,'tickdir','out')
