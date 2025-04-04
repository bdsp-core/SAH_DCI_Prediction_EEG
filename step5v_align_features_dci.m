% viasualization: heatmap and difference analysis of two groups

clear
close all;

% for vascular and asymmetry:
% load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_vas');
% features = vas_features_dci;
% feature_score = vas_feature_score_dci;
% features = asym_features_dci;
% feature_score = asym_feature_score_dci;
% ifeature = 2;
% ivascular = 1;
% feature_all = features{ifeature,ivascular};
% score_all = feature_score{ifeature,ivascular};

% for global features:
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_global.mat');
features = features_dci;%spike_dci
feature_score = feature_score_dci; % spike_score
ifeature = 2;
feature_all = features{ifeature};
score_all = feature_score{ifeature};

load('no_dci_match.mat');
pts_info = readtable('sah_dci_FINALadjud.csv');
time_shift_all = zeros(size(dci_labels));
pts_info.SID_str = arrayfun(@(x) sprintf('%03d', x), pts_info.SID, 'UniformOutput', false);


% get time info from chart review
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
        time_shift_all(ipatient) = h; % time after bleed to dci onset
    end
end
% match time shift for non dci patient with time shift of dci patient in
% the same position
for i = 1:length(no_dci_index)
    time_shift_all(no_dci_index(i)) = time_shift_all(no_dci_match(i));
end

%% Visualization

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
vascular_save = {'ACA_left'
    'MCA_left'
    'PCA_left'
    'ACA_right'
    'MCA_right'
    'PCA_right'
    };

% feature_all = 10*log10(feature_all);
index_tmp = find(dci_labels==0);
feature_good = feature_all(dci_labels==0,:);
feature_mean = nanmean(feature_good');
[a,index_good]=sort(feature_mean); % sorts them in ascending order.
index_good = index_tmp(index_good);
feature_good = feature_all(index_good,:);
unique_names_good = unique_names_dci(index_good);
score_good = score_all(index_good,:);

index_tmp = find(dci_labels==1);
feature_bad = feature_all(dci_labels==1,:);
feature_mean = nanmean(feature_bad');
[a,index_bad]=sort(feature_mean); % sorts them in ascending order.
index_bad = index_tmp(index_bad);
feature_bad = feature_all(index_bad,:);
unique_names_bad = unique_names_dci(index_bad);
score_bad = score_all(index_bad,:);

unique_names1 = [unique_names_good';unique_names_bad'];
feature_all1 = [feature_good;feature_bad];
score_all = [score_good;score_bad];


fig = figure('units','normalized','outerposition',[0 0 1 1]);
mn = prctile(feature_all1(:),1); 
mx = prctile(feature_all1(:),99);% 90
imagesc(feature_all1(:,700:1300),[mn,mx]);
% Only shows the time range 700:1300 (centered around DCI occurrence).
colormap;
set(gca, 'YTick', (1/num_pts:1/num_pts:1)*num_pts, 'YTickLabel', unique_names1)
set(gca, 'XTick', 1:50:601, 'XTickLabel', cellstr(num2str(((700:50:1300)-dci_point)')));
xlabel('Hours After DCI');
% title([feature_names{ifeature},'\_',vascular_names{ivascular}]);
title(feature_names{ifeature});
% title('Spike Frequency')
colorbar;

%% Mean and Standard Deviation Plot
figure
mean_feature = nanmean(feature_good);
std_feature = nanstd(feature_good);
valid_pts = sum(~isnan(feature_good));
range_start = 725;
range_end = 1200-7;
xb= range_start:range_end;
yb = mean_feature(range_start:range_end);
[l,p] = boundedline(range_start:range_end, mean_feature(range_start:range_end), std_feature(range_start:range_end)./sqrt(valid_pts(range_start:range_end)),'--b*','alpha');% --b* --ro
h1 = outlinebounds(l,p);

mean_feature = nanmean(feature_bad);
std_feature = nanstd(feature_bad);
valid_pts = sum(~isnan(feature_bad));
xr = range_start:range_end;
yr = mean_feature(range_start:range_end);
[l,p] = boundedline(range_start:range_end, mean_feature(range_start:range_end), std_feature(range_start:range_end)./sqrt(valid_pts(range_start:range_end)),'--ro','alpha');% --b* --ro
h2 = outlinebounds(l,p);

xlim([range_start-1,range_end]),xlabel('Hours After DCI')
set(gca, 'XTick', range_start+25:50:range_end, 'XTickLabel', cellstr(num2str(((range_start+25:50:range_end)-dci_point)')));
ylabel('Feature Values')
title(feature_names{ifeature});
% title([feature_names{ifeature},'\_',vascular_names{ivascular}]);
% title('Spike Frequency');
axis tight;
legend([h1 h2],{'No DCI','DCI'});
% saveas(gcf,['D:\Research\SAH_DCI_Jenn\results\alignWithDCI\vascular_territory\', feature_names{ifeature},'_',vascular_save{ivascular},'.png']);
% close all;

