% This file is for spike processing

clear;
file_folder = 'Z:\Projects\Weilong\SAH_DCI\Spikes_New_SSD\';
file_folder_mask = 'Z:\Projects\Weilong\SAH_DCI\Preprocessed\';
save_folder = 'F:\SAH_DCI\Spikes_New_SSD\PerHour30\';
pts_info = readtable('D:\Research\SAH_DCI_Jenn\sah_dci_FINALadjud.csv');

list = dir([file_folder_mask,'*.mat']);
isfile=~[list.isdir]; %determine index of files vs folders
filenames={list(isfile).name}; %create cell array of file names
num_file = length(filenames);
num_channel = 18;
win_feature = 5; % min
spike_threshold = 0.30; %0.43 0.5
for i = 1:num_file
    underlineLocations = find(filenames{i} == '_');
    file_name{i} = filenames{i}(underlineLocations(2)-3:underlineLocations(2)-1);
    file_date{i} = filenames{i}(underlineLocations(end-1)+1:underlineLocations(end)-1);
    file_time{i} = filenames{i}(underlineLocations(end)+1:underlineLocations(end)+6);
    dateRecord = datetime(str2num(file_date{i}(1:4)),str2num(file_date{i}(5:6)),str2num(file_date{i}(7:8)),str2num(file_time{i}(1:2)),str2num(file_time{i}(3:4)),str2num(file_time{i}(5:6)));
    
    indx = ~cellfun(@isempty, strfind(pts_info.SID,file_name{i}));
    indx_tmp = find(indx==true);
    dates = datetime(char(pts_info.DayOfBleed(indx_tmp(1))),'Format','MM/dd/yyyy');
    times = datenum(char(pts_info.TimeOfBleed(indx_tmp(1))));
    dates = datetime(dates,'Format','MM/dd/yyyy HH:mm:SS');
    times = datetime(times,'ConvertFrom','datenum','Format','MM/dd/yyyy HH:mm:SS');
    timeBleed = dates+timeofday(times);
    dt_start = between(timeBleed,dateRecord,'time');
    
%     try
    load([file_folder_mask,filenames{i}],'seg_masks');
    load([file_folder,filenames{i}],'yp','artifact');
    file_length = size(seg_masks,1);
    durationRecord = seconds(file_length*5);
    dt_end = dt_start+durationRecord;
    
    dt_start_vec = datevec(dt_start);
    dt_end_vec = datevec(dt_end);
    dhour_start = dt_start_vec(4)+1;
    dhour_end = dt_end_vec(4);
    num_hour = dhour_end-dhour_start+1;
    
    fs_spike = 128;
    yp(artifact==1)=NaN;
    len = floor(length(yp)/(5*fs_spike));
    yp = reshape(yp(1:len*(5*fs_spike)),[(5*fs_spike),len]);
    
%     eeg_pos = {};
    eeg_masks = {};
    eeg_spike = {};
%     eeg_specs = {};
%     eeg_data = {};
    eeg_hour = {};
%     eeg_feature = {};
                
    if num_hour>0
        for j = 1:num_hour
            current_timediff = hours(dt_start_vec(4)+j)-dt_start; 
            if j==1
                current_duration = current_timediff;
            else
                current_duration = hours(1);
            end
            
            current_timediff = datevec(current_timediff);
            current_duration = datevec(current_duration);
            
            current_indexdiff = round((((current_timediff(3)*24+current_timediff(4))*60+current_timediff(5))*60+current_timediff(6))/5);
            current_index_duration = round((((current_duration(3)*24+current_duration(4))*60+current_duration(5))*60+current_duration(6))/5);
            
            num_win = floor(current_index_duration*5/(win_feature*60));
            %current_eeg = EEG_segs(current_indexdiff-num_win*win_feature*60/5+1:current_indexdiff,:,:);
            %current_masks = seg_masks(current_indexdiff-num_win*win_feature*60/5+1:current_indexdiff,:);
            
            for sliding_win = 1:num_win
                time_range = current_indexdiff-sliding_win*win_feature*60/5+1:current_indexdiff-(sliding_win-1)*win_feature*60/5;
                %artifact_check = sum(cellfun(@isempty, strfind(string(seg_masks(time_range,:)),'normal')));

                eeg_masks{j,sliding_win} = seg_masks(time_range,:);
                eeg_hour{j,sliding_win} = dhour_start+j-1;
                
                spike_segs = yp(:,time_range);
                spike_segs = reshape(spike_segs,[size(spike_segs,1)*size(spike_segs,2),1]);
                spike_timepoints = find(spike_segs>=spike_threshold);
                spike_ranges = convert_indices_to_index_ranges(spike_timepoints);
                eeg_spike{j,sliding_win} = size(spike_ranges,1);
            end
            
        end
        save([save_folder,filenames{i}],'eeg_masks','eeg_spike','eeg_hour');
    else
        current_timediff = dt_end-dt_start; 
            current_duration = current_timediff;
            
            current_timediff = datevec(current_timediff);
            current_duration = datevec(current_duration);
            
            current_indexdiff = round((((current_timediff(3)*24+current_timediff(4))*60+current_timediff(5))*60+current_timediff(6))/5);
            current_index_duration = round((((current_duration(3)*24+current_duration(4))*60+current_duration(5))*60+current_duration(6))/5);
            
            num_win = floor(current_index_duration*5/(win_feature*60));
            %current_eeg = EEG_segs(current_indexdiff-num_win*win_feature*60/5+1:current_indexdiff,:,:);
            %current_masks = seg_masks(current_indexdiff-num_win*win_feature*60/5+1:current_indexdiff,:);
            j = 1;
            eeg_masks = cell(1,12);
            eeg_hour = cell(1,12);
            eeg_spike = cell(1,12);
            
            for sliding_win = 1:num_win
                time_range = current_indexdiff-sliding_win*win_feature*60/5+1:current_indexdiff-(sliding_win-1)*win_feature*60/5;
                %artifact_check = sum(cellfun(@isempty, strfind(string(seg_masks(time_range,:)),'normal')));

                eeg_masks{j,sliding_win} = seg_masks(time_range,:);
                eeg_hour{j,sliding_win} = dhour_start+j-1;
                
                spike_segs = yp(:,time_range);
                spike_segs = reshape(spike_segs,[size(spike_segs,1)*size(spike_segs,2),1]);
                spike_timepoints = find(spike_segs>=spike_threshold);
                spike_ranges = convert_indices_to_index_ranges(spike_timepoints);
                eeg_spike{j,sliding_win} = size(spike_ranges,1);

            end
            save([save_folder,filenames{i}],'eeg_masks','eeg_spike','eeg_hour');
    end

end