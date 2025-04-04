% This file extract the feature from step 2 to 7 global features:
% 'Shannon', 'AlphaDeltaRatio', 'TotalPow', 'DeltaPow', 'ThetaPow', 'AlphaPow', 'PerAlphaVar'
% calculate feature_score - the weights based on artifact detection
% and save the features into one single file 'all_features_NEWspike30.mat'

% Folder 'Features' should have all patients' data, in this example
% 'all_features_NEWspike30.mat' is already generated

clear
file_folder = 'E:\Zhongwei\SAH Code Publish\Features\';
% spike_folder = 'F:\SAH_DCI\Spikes_New_SSD\PerHour30\';
save_folder = 'E:\Zhongwei\SAH Code Publish\Extract_Features\';
pts_info = readtable('E:\Zhongwei\SAH Code Publish\sah_dci_FINALadjud.csv');

list = dir([file_folder,'*.mat']);
isfile=~[list.isdir]; %determine index of files vs folders
filenames={list(isfile).name}; %create cell array of file names
num_file = length(filenames);

%% find longest recording among all patients to figure out max matrix size
max_hour = 0;
for i = 1:num_file
    underlineLocations = find(filenames{i} == '_');
    pt_name{i} = filenames{i}(underlineLocations(2)-3:underlineLocations(2)-1);
    load ([file_folder,filenames{i}],'eeg_hour');
    if max_hour<eeg_hour{end,1}
        max_hour = eeg_hour{end,1};
        max_filename = filenames{i};
    end
end
%%

feature_names = {'Shannon'
    'AlphaDeltaRatio'
    'TotalPow'
    'DeltaPow'
    'ThetaPow'
    'AlphaPow'
    'PerAlphaVar'
    };
[unique_names, idx ,idx2] = uniquecell(pt_name);
num_pts = length(unique_names);
feature_all = NaN(num_pts,max_hour);
score_all = NaN(num_pts,max_hour);
dci_labels = zeros(num_pts,1);


pts_info.SID_str = arrayfun(@(x) sprintf('%03d', x), pts_info.SID, 'UniformOutput', false);
for feature_num = 1:length(feature_names) % change this number to calculate each feature as listed
for i = 1:num_pts
    current_pts = unique_names{i};
    indx = contains(pts_info.SID_str,current_pts);
    indx_tmp = find(indx==true);
    dci_date = pts_info.Final_DCIDate(indx_tmp(1));

    if contains(string(dci_date), 'NaN')
        dci_labels(i) = 0;
    else
        dci_labels(i) = 1;
    end
    unique_names_dci{i} = [unique_names{i},'(',num2str(dci_labels(i)),')'];
    file_index = find(idx2==i);
    
    for j = 1:length(file_index)
        
        load([file_folder,filenames{file_index(j)}],'eeg_hour','eeg_feature','eeg_masks');
        for k = 1:size(eeg_hour,1)
            index = cellfun(@(x)isequal(x,eeg_hour{k,1}),eeg_hour);
            [row,col] = find(index);
            
            current_mask = eeg_masks(k,col);
            current_feature = eeg_feature(k,col);
            weight_clean = zeros(1,length(col));
            feature_raw = zeros(1,length(col));
            
            %% calculating quality scores for artifacts
            for m = 1:length(col)
                indx = ~cellfun(@isempty, strfind(string(current_mask{m}),'normal'));
                weight_clean(m) = sum(indx)/length(indx);
                
                feature_raw(m) = mean(rmoutliers(current_feature{m}(feature_num:7:end)));
                
            end
            weight_clean_norm = weight_clean/sum(weight_clean);
            feature_clean = sum(weight_clean_norm.*feature_raw);
            
            feature_all(i,eeg_hour{k,1}) = feature_clean;
            score_all(i,eeg_hour{k,1}) = max(weight_clean);
        end
    end
end
    features{feature_num} = feature_all;
    feature_score{feature_num} = score_all;
end


%% spike
feature_all = NaN(num_pts,max_hour);
score_all = NaN(num_pts,max_hour);
for i = 1:num_pts
    current_pts = unique_names{i};
    indx = contains(pts_info.SID_str,current_pts);
    indx_tmp = find(indx==true);
    dci_date = pts_info.Final_DCIDate(indx_tmp(1));
    if contains(string(dci_date),'NaN')
        dci_labels(i) = 0;
    else
        dci_labels(i) = 1;
    end
    unique_names_dci{i} = [unique_names{i},'(',num2str(dci_labels(i)),')'];
    file_index = find(idx2==i);
    
    for j = 1:length(file_index)
        
        load([spike_folder,filenames{file_index(j)}],'eeg_hour','eeg_spike','eeg_masks');
        for k = 1:size(eeg_hour,1)
            index = cellfun(@(x)isequal(x,eeg_hour{k,1}),eeg_hour);
            [row,col] = find(index);
            
            current_mask = eeg_masks(k,col);
            current_feature = eeg_spike(k,col);
            weight_clean = zeros(1,length(col));
            feature_raw = zeros(1,length(col));
            
            for m = 1:length(col)
                indx = contains(string(current_mask{m}),'normal');
                weight_clean(m) = sum(indx)/length(indx);
                feature_raw(m) = mean(mean(current_feature{m}));
            end
            weight_clean_norm = weight_clean/sum(weight_clean);
            feature_clean = sum(weight_clean_norm.*feature_raw);%feature_clean = sum(feature_raw);%
            
            feature_all(i,eeg_hour{k,1}) = feature_clean;
            score_all(i,eeg_hour{k,1}) = max(weight_clean);
        end
    end
end
spike = feature_all;
spike_score = score_all;

save([save_folder,'all_features_NEWspike30'],'features','feature_score','spike','spike_score','dci_labels','unique_names','unique_names_dci');

