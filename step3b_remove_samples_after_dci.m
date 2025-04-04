% set all values after DCI onset to 'NaN

clear
close all;
load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\all_features_NEWspike30.mat');

pts_info = readtable('sah_dci_FINALadjud.csv');
dci_time_after_bleed = zeros(1,length(unique_names));

% Convert numeric SID to strings with zero-padding (assuming 3-digit format)
pts_info.SID_str = arrayfun(@(x) sprintf('%03d', x), pts_info.SID, 'UniformOutput', false);

%%
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
        % hms extracts the hours, minutes, and seconds from the duration object and returns them as three separate outputs.
        dci_time_after_bleed(ipatient) = h;
    end
end

for i = 1:size(spike,1)
    if ~isnan(dci_time_after_bleed(i)) && dci_time_after_bleed(i) ~= 0
        spike(i,dci_time_after_bleed(i):end) = NaN;
        spike_score(i,dci_time_after_bleed(i):end) = NaN;
        for j = 1:7
            features{j}(i,dci_time_after_bleed(i):end) = NaN;
            feature_score{j}(i,dci_time_after_bleed(i):end) = NaN;
        end
    end
end
save('.\ExtractFeatures\features_remove_dci_global.mat','dci_labels','feature_score','features','spike','spike_score','unique_names','unique_names_dci');

%% asymmetry

load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\asymmetry_subtraction_features.mat');

for i = 1:length(unique_names)
    if ~isnan(dci_time_after_bleed(i)) && dci_time_after_bleed(i) ~= 0
        for ifea = 1:size(asym_features,1)
            for jfea = 1:size(asym_features,2)
                asym_features{ifea,jfea}(i,dci_time_after_bleed(i):end) = NaN;
            end
        end
    end
end

save('.\ExtractFeatures\features_remove_dci_asym.mat','asym_features', 'dci_labels', 'feature_score', 'unique_names','unique_names_dci');



%% vascular

load('E:\Zhongwei\SAH Code Publish\ExtractFeatures\vascular_features.mat');

for i = 1:length(unique_names)
    if ~isnan(dci_time_after_bleed(i)) && dci_time_after_bleed(i) ~= 0
        for ifea = 1:size(features,1)
            for jfea = 1:size(features,2)
                features{ifea,jfea}(i,dci_time_after_bleed(i):end) = NaN;
                feature_score{ifea,jfea}(i,dci_time_after_bleed(i):end) = NaN;
            end
        end
    end
end

save('E:\Zhongwei\SAH Code Publish\ExtractFeatures\features_remove_dci_vas.mat','dci_labels','features','feature_score','features','unique_names','unique_names_dci');