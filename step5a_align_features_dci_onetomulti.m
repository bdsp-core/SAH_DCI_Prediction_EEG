% match one DCI patient to multiple non-DCI patient(jenn's code; method of paper)
clear
close all;

file_vas = 'E:\Zhongwei\SAH Code Publish\ExtractFeatures\vascular_features.mat';
file_asym = 'E:\Zhongwei\SAH Code Publish\ExtractFeatures\asymmetry_subtraction_features.mat';
file_glob = 'E:\Zhongwei\SAH Code Publish\ExtractFeatures\all_features_NEWspike30.mat';

% Things to do:
file = file_glob;

load(file);
load('ctrl_matchtime_dci.mat');
total_cases_after_dci = length(cat(1,match_ctrls{:}))+length(ctrl_sids);
%% 
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

%%
if strcmp(file, file_asym)
    features = asym_features;
end

total_hour = size(features{1,1},2);
dci_point = 1000;
for i = 1:size(features,1)
    for j = 1:size(features,2)
        current_feature = features{i,j};
        fea_temp = NaN(size(pts_index,2),2000);
        for ipatient = 1:length(pts_index)
            shift_length = dci_point-time_shift_all(ipatient);
            fea_temp(ipatient,(1:total_hour)+shift_length) = current_feature(pts_index(ipatient),:);
        end
    features_dci{i,j} = fea_temp;
    end
end

if strcmp(file, file_glob)
    current_feature = spike;
    fea_temp = NaN(size(pts_index,2),2000);
    for ipatient = 1:length(pts_index)
        shift_length = dci_point-time_shift_all(ipatient);
        fea_temp(ipatient,(1:total_hour)+shift_length) = current_feature(pts_index(ipatient),:);
    end
    spike_dci = fea_temp;

    save('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_global','features_dci','spike_dci','dci_labels_new','unique_names','unique_names_dci');
elseif strcmp(file, file_asym)
    asym_features_dci = features_dci;
    save('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_asym','asym_features_dci', 'dci_labels_new','unique_names','unique_names_dci');
else
    vas_features_dci = features_dci;    
    save('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_align_dci_jenn_vas','vas_features_dci','dci_labels_new','unique_names','unique_names_dci');
end


