%% 
% Step1: convert Artifact reduced edf file to a mat file, 
% Read and resample ARed file to 200, 
% Save available channels and segment the EEG siganls after notch and bandpass filter
% 
% This MATLAB script processes preprocessed EEG data files of step1,
% extracts features, and saves the results in a structured format. 


clear;

% Things to Do
% folder path
file_folder = 'E:\Zhongwei\SAH Code Publish\Preprocessed\';
save_folder = 'E:\Zhongwei\SAH Code Publish\Features\';

% csv where pts time of bleed info saved
pts_info = readtable('E:\Zhongwei\SAH Code Publish\sah_dci_FINALadjud.csv');


pts_info.SID_str = arrayfun(@(x) sprintf('%03d', x), pts_info.SID, 'UniformOutput', false);

list = dir([file_folder,'*.mat']);
isfile=~[list.isdir]; %determine index of files vs folders
filenames={list(isfile).name}; %create cell array of file names
% filenames = flip(filenames,2);
num_file = length(filenames);
num_channel = 18;
win_feature = 5; % min
for i = 1:num_file

    % Extract the time of Bleed in format mm/dd/yyyy hh:mm:ss using the
    % record in csv under the name of current .mat file. Calculate the
    % time elapse before bleeding: dt_start = [recording start:time of bleed]
    underlineLocations = find(filenames{i} == '_');
    file_name{i} = filenames{i}(underlineLocations(2)-3:underlineLocations(2)-1);
    file_date{i} = filenames{i}(underlineLocations(end-1)-8:underlineLocations(end)-1);
    file_time{i} = filenames{i}(underlineLocations(end)-6:underlineLocations(end)-1);
    dateRecord = datetime(str2num(file_date{i}(1:4)),str2num(file_date{i}(5:6)),str2num(file_date{i}(7:8)),str2num(file_time{i}(1:2)),str2num(file_time{i}(3:4)),str2num(file_time{i}(5:6)));
    
    indx = contains(pts_info.SID_str,file_name{i});
    indx_tmp = find(indx==true);
    dates = datetime(char(pts_info.DayOfBleed(indx_tmp(1))),'Format','MM/dd/yyyy');
    times = datenum(char(pts_info.TimeOfBleed(indx_tmp(1))));
    dates = datetime(dates,'Format','MM/dd/yyyy HH:mm:SS');
    times = datetime(times,'ConvertFrom','datenum','Format','MM/dd/yyyy HH:mm:SS');
    timeBleed = dates+timeofday(times);
    dt_start = between(timeBleed,dateRecord,'time');

    load([file_folder,filenames{i}],'EEG_segs_bipolar','EEG_specs','EEG_frequency','seg_masks','Fs','channel_names','burst_suppression');
    EEG_segs = EEG_segs_bipolar; % (x,18,1000) = (windows_num, channel_num, window_size); windowsize = number of samples in each window (window_time * Fs)
    clear EEG_segs_bipolar
    file_length = size(EEG_segs,1);
    durationRecord = seconds(file_length*5);

    dt_end = dt_start+durationRecord;

    dt_start_vec = datevec(dt_start);
    dt_end_vec = datevec(dt_end);
    dhour_start = dt_start_vec(4)+1;
    dhour_end = dt_end_vec(4);
    num_hour = dhour_end-dhour_start+1;


    % The use of dt_end is not clear, seems like num_hour is rounded equal
    % to durationRecord

    % Segments EEG into 5-minute windows.
    % Extracts artifact masks, burst suppression data, and EEG power spectra.
    % Calls a function Feature_Subset to extract features from each EEG segment.
    eeg_pos = {};
    eeg_masks = {};
    eeg_bs = {};
    eeg_specs = {};
    eeg_data = {};
    eeg_hour = {};
    eeg_feature = {};
                
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
            
            for sliding_win = 1:num_win
                time_range = current_indexdiff-sliding_win*win_feature*60/5+1:current_indexdiff-(sliding_win-1)*win_feature*60/5;

                eeg_pos{j,sliding_win} = [time_range(1) time_range(end)];
                eeg_masks{j,sliding_win} = seg_masks(time_range,:);
                eeg_bs{j,sliding_win} = burst_suppression(time_range,:);
                eeg_specs{j,sliding_win} = EEG_specs(time_range,:);
                eeg_data{j,sliding_win} = EEG_segs(time_range,:,:);
                eeg_hour{j,sliding_win} = dhour_start+j-1;
                
                eeg = EEG_segs(time_range,:,:);
                eeg = permute(eeg,[3,1,2]);
                eeg = reshape(eeg,[size(eeg,1)*size(eeg,2),size(eeg,3)]);
                [num_sample, num_chan] = size(eeg);
                feature =[];
                for ichannel = 1:num_chan
                    [out, header] = Feature_Subset( eeg, ichannel, Fs );
                    feature = [feature out];
                end
                
                eeg_feature{j,sliding_win} = feature;

            end
            
        end
        save([save_folder,filenames{i}],'eeg_pos','eeg_masks','eeg_data','eeg_bs','eeg_specs','eeg_feature','Fs','channel_names','EEG_frequency','eeg_hour');
    else
        current_timediff = dt_end-dt_start; 
            current_duration = current_timediff;
            
            current_timediff = datevec(current_timediff);
            current_duration = datevec(current_duration);
            
            current_indexdiff = round((((current_timediff(3)*24+current_timediff(4))*60+current_timediff(5))*60+current_timediff(6))/5);
            current_index_duration = round((((current_duration(3)*24+current_duration(4))*60+current_duration(5))*60+current_duration(6))/5);
            
            num_win = floor(current_index_duration*5/(win_feature*60));

            j = 1;
            eeg_pos = cell(1,12);
            eeg_masks = cell(1,12);
            eeg_bs = cell(1,12);
            eeg_specs = cell(1,12);
            eeg_data = cell(1,12);
            eeg_hour = cell(1,12);
            eeg_feature = cell(1,12);
            
            for sliding_win = 1:num_win
                time_range = current_indexdiff-sliding_win*win_feature*60/5+1:current_indexdiff-(sliding_win-1)*win_feature*60/5;

                eeg_pos{j,sliding_win} = [time_range(1) time_range(end)];
                eeg_masks{j,sliding_win} = seg_masks(time_range,:);
                eeg_bs{j,sliding_win} = burst_suppression(time_range,:);
                eeg_specs{j,sliding_win} = EEG_specs(time_range,:);
                eeg_data{j,sliding_win} = EEG_segs(time_range,:,:);
                eeg_hour{j,sliding_win} = dhour_start-1;
                
                eeg = EEG_segs(time_range,:,:);
                eeg = permute(eeg,[3,1,2]);
                eeg = reshape(eeg,[size(eeg,1)*size(eeg,2),size(eeg,3)]);
                [num_sample, num_chan] = size(eeg);
           
                feature =[];
                for ichannel = 1:num_chan
                    [out, header] = Feature_Subset( eeg, ichannel, Fs );
                    feature = [feature out];
                end
                
                eeg_feature{j,sliding_win} = feature;

            end
        save([save_folder,filenames{i}],'eeg_pos','eeg_masks','eeg_data','eeg_bs','eeg_specs','eeg_feature','Fs','channel_names','EEG_frequency','eeg_hour');
    end
end