% This file is for spike processing

clear
load('compiled_spks_persystpercept_EEGtimes30d_v3.mat')
max_hour = 584;
spike = nan(length(compiled_spks_persystpercept),max_hour); % maximum hour 600
hour_index = (1:max_hour)*3600;
for i = 1:length(compiled_spks_persystpercept)
    spike_pro = compiled_spks_persystpercept{i};
    for j = 1:max_hour
        index = find(spike_pro(:,1)>(j-1)*3600+1&spike_pro(:,1)<=j*3600);
        if ~isempty(index)
            spike(i,j) = length(spike_pro(index,2));
        end
    end
end