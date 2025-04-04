% generate global/vascular/asymmetry features aligned with DCI onset.

clear
close all;

file_vas = 'E:\Zhongwei\SAH Code Publish\ExtractFeatures\vascular_features.mat';
file_asym = 'E:\Zhongwei\SAH Code Publish\ExtractFeatures\asymmetry_subtraction_features.mat';
file_glob = 'E:\Zhongwei\SAH Code Publish\ExtractFeatures\all_features_NEWspike30.mat';

% Things to do:
file = file_asym;

load(file);
load('no_dci_match.mat');
pts_info = readtable('sah_dci_FINALadjud.csv');
time_shift_all = zeros(size(dci_labels));
pts_info.SID_str = arrayfun(@(x) sprintf('%03d', x), pts_info.SID, 'UniformOutput', false);


%% get time info from chart review
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

%% align global/asym/vas features
if strcmp(file, file_asym)
    features = asym_features;
end

dci_point = 1000; % each pt will be aligned to dci onset time point=1000
total_hour = size(features{1,1},2); % total_hour=584
for i = 1:size(features,1) % 1:7
    for j = 1:size(features,2) % 1:6
    current_feature = features{i,j};
    current_scores = feature_score{i,j};
    fea_temp = NaN(size(features{1,1},1),2000); % 113*2000
    score_temp = NaN(size(features{1,1},1),2000);
    for ipatient = 1:length(unique_names)
        shift_length = dci_point-time_shift_all(ipatient);
        % each pt will be aligned to dci onset time point=1000
        % for non dci pt, shift length will be the same as dci pt
        fea_temp(ipatient,(1:total_hour)+shift_length) = current_feature(ipatient,:);
        score_temp(ipatient,(1:total_hour)+shift_length) = current_scores(ipatient,:);
    end
    features_dci{i,j} = fea_temp;
    feature_score_dci{i,j} = score_temp;
    end
end

if strcmp(file, file_glob)
    current_feature = spike;
    current_scores = spike_score;
    fea_temp = NaN(size(spike,1),2000);
    score_temp = NaN(size(spike,1),2000);
    for ipatient = 1:length(unique_names)
        shift_length = dci_point-time_shift_all(ipatient);
        fea_temp(ipatient,(1:total_hour)+shift_length) = current_feature(ipatient,:);
        score_temp(ipatient,(1:total_hour)+shift_length) = current_scores(ipatient,:);
    end
    spike_dci = fea_temp;
    spike_score = score_temp;

    save('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_global','features_dci','feature_score_dci','spike_dci','spike_score','dci_labels','unique_names','unique_names_dci');
elseif strcmp(file, file_asym)
    asym_features_dci = features_dci;
    asym_feature_score_dci = feature_score_dci;

    % save('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_asym_rmDCI','asym_features_dci','asym_feature_score_dci','dci_labels','unique_names','unique_names_dci');
    save('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_asym','asym_features_dci','asym_feature_score_dci','dci_labels','unique_names','unique_names_dci');
else
    vas_features_dci = features_dci;
    vas_feature_score_dci = feature_score_dci;
    
    save('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_vas','vas_features_dci','vas_feature_score_dci','dci_labels','unique_names','unique_names_dci');
end

