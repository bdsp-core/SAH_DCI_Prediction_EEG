% Generate vascular features from edf files
% Generate asymmetry features based on vascular features
clear

load('utile.mat')
%%
file_folder = 'E:\Zhongwei\SAH Code Publish\Features';
save_folder = 'E:\Zhongwei\SAH Code Publish\ExtractFeatures';
pts_info = readtable('sah_dci_FINALadjud_outcomes_withclinic_final.csv');

%% 
vascular_names = {'ACA_left'
    'MCA_left'
    'PCA_left'
    'ACA_right'
    'MCA_right'
    'PCA_right'
    };
vascular_index = {[1,9],...
    [10,11,2,3],...
    [12,4],...
    [13,5],...
    [14,15,6,7],...
    [16,8]
    };

num_pts = numel(unique_names);
feature_all = NaN(num_pts,max_hour);
score_all = NaN(num_pts,max_hour);

for feature_num = 1:7
    for j = 1:numel(vascular_names)
        features{feature_num,j} = NaN(num_pts,max_hour);
        feature_score{feature_num,j} = NaN(num_pts,max_hour);
    end
end


%% calculate & merge features
% Initialize waitbar
h = waitbar(0, 'Merging features...');

% Process each feature
for i = 1:num_pts
% for i = 16:18


    % Get current patient ID and find matching rows in pts_info
    waitbar(i / num_pts, h, sprintf('Processing patient %d of %d...', i, num_pts));
    current_pt = strip(unique_names{i}, 'left', '0')
    
    % Find and process files associated with the current patient
    file_idx = find(idx2 == i);
    for j = 1:numel(file_idx)
        data = load(fullfile(file_folder, filenames{file_idx(j)}), 'eeg_hour', 'eeg_feature', 'eeg_masks');
        % Skip files with missing data fields
        if ~all(isfield(data, {'eeg_hour', 'eeg_feature', 'eeg_masks'}))
            warning('Missing expected data in file %s.', filenames{file_idx(j)});
            continue;
        end
        
        for k = 1:size(data.eeg_hour,1)
%             k = k_tmp;
%             if data.eeg_hour{k_tmp,1}==[]
%                 k = k-1;
%                 continue;
%             end
            index = cellfun(@(x)isequal(x,data.eeg_hour{k,1}),data.eeg_hour);
            [row,col] = find(index);
            
            current_mask = data.eeg_masks(k,col);
            current_feature = data.eeg_feature(k,col);
            weight_clean = zeros(1,length(col));
            feature_raw = zeros(1,length(col));

            for feature_num = 1:numel(feature_names)
                for ivascular = 1:numel(vascular_names)
            % calculating quality scores for artifacts
                
                    for m = 1:length(col)
                        indx = contains(string(current_mask{m}),'normal');
                        weight_clean(m) = sum(indx)/length(indx);
                        
                        % feature_raw is the mean value of one feature(out of
                        % 7) at that col
                        current_fea = current_feature{m}(feature_num:7:end);
                        if ~isempty(current_fea)
                            feature_raw(m) = mean(current_fea(vascular_index{ivascular}));
                        end
                    end
    
                    weight_clean_norm = weight_clean/sum(weight_clean);
                    feature_clean = sum(weight_clean_norm.*feature_raw);
    
                    features{feature_num,ivascular}(i,data.eeg_hour{k,1}) = feature_clean;
                    feature_score{feature_num,ivascular}(i,data.eeg_hour{k,1}) = max(weight_clean);
                end
            end
        end
    end
end

close(h);

save(fullfile(save_folder,'vascular_features'),'features','feature_score','dci_labels','unique_names','unique_names_dci', '-v7.3');

%% generate asymmetry feature

load(fullfile(save_folder,'vascular_features.mat'))

% asym_features is the ACA, MCA, PCA absolute difference of left vs right
asym_features = cell(7,3);
for ifeature = 1:7
    for ivascular = 1:3
        feature_all = abs(features{ifeature,ivascular}-features{ifeature,ivascular+3});
        % feature_all = abs((10*log10(features{ifeature,ivascular}))-(10*log10(features{ifeature,ivascular+3})));
        % feature_all = max(features{ifeature,ivascular}./features{ifeature,ivascular+3},features{ifeature,ivascular+3}./features{ifeature,ivascular});
        % feature_all = max((10*log10(features{ifeature,ivascular}))./(10*log10(features{ifeature,ivascular+3})),(10*log10(features{ifeature,ivascular+3}))./(10*log10(features{ifeature,ivascular})));
        % feature_all = features{ifeature,ivascular}./features{ifeature,ivascular+3};

        asym_features{ifeature,ivascular} = feature_all;

    end
end

save(fullfile(save_folder,'asymmetry_subtraction_features'),'asym_features','feature_score','dci_labels','unique_names','unique_names_dci');





%% old method for pav 6h calculation, for verification, not used
%% ATR to hourly

atr_all = NaN(12, max_hour); % max windows per hour is 12 (for 5 mins)

feature_num = 7;

for ivascular = 1:length(vascular_names)
    for i = 1

        file_index = find(idx2==i);
        
        for j = 1:length(file_index)
            load(fullfile(file_folder,filenames{file_index(j)}),'eeg_hour','eeg_feature','eeg_masks');
            
            if isempty(eeg_hour{1,1})
                disp(['first eeg_hour empty in ', filenames{file_index(j)}]);
                index = cellfun(@(x)isequal(x,eeg_hour{1,1}),eeg_hour)
                [row,col] = find(index)
%                 for tmp = 1:size(eeg_hour, 2) % Loop through each column in the first row
%                     disp(isempty(eeg_hour{1, tmp}));
%                 end
                eeg_hour = eeg_hour(2:end,:);
            end
            for k = 1:size(eeg_hour,1)
                if isempty(eeg_hour{k,1})
                    disp(['eeg_hour empty', filenames{file_index(j)},k]);
                    pause;
                end

                index = cellfun(@(x)isequal(x,eeg_hour{k,1}),eeg_hour);
                [row,col] = find(index);
                
                current_feature = eeg_feature(k,col);
                feature_raw = zeros(1,length(col));
                
                %% calculating quality scores for artifacts
                for m = 1:length(col)
                    current_fea = current_feature{m}(feature_num:7:end);
                    if ~isempty(current_fea)
                        feature_raw(m) = mean(current_fea(vascular_index{ivascular}));
                    end
                    
                end
                feature_clean = feature_raw;%feature_clean = sum(feature_raw);%
                
                if ~(sum(isnan(feature_clean)) == length(col))
                    atr_all(1:length(feature_clean), eeg_hour{k,1}) = feature_clean';
                end
                
            end
        end
        
        %% now calculate PAV for every 6h
        
        % the first PAV value will be at hour 6 (taking into account hours 1-6)
        pav = NaN(1,max_hour);
        for hr = 6:max_hour
            if ~isnan(mean(atr_all(:, hr)))
                pav(hr) = iqr(atr_all(:, hr-5:hr), 'all');
            end
        end
        
        feature_all(i, :) = pav;
    end
    
    features{feature_num,ivascular} = feature_all;
end

