#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

This file is modified by Zhongwei. Unused code is deleted for readability,
refer og code if encountered super suspicious problem.
"""

from collections import Counter
import numpy as np
from joblib import Parallel, delayed
from scipy.signal import detrend
import matplotlib.pyplot as plt
import mne
mne.set_log_level(verbose='WARNING')
from mne.filter import filter_data, notch_filter



def segment_EEG_without_detection(EEG, Ch_names, window_time, step_time, Fs, notch_freq=None, bandpass_freq=None, start_end_remove_window_num=0,
                amplitude_thres=500, n_jobs=1, to_remove_mean=False):#
    """Segment EEG signals, similar to segment_EEG.py, but does not apply artifact detection.
    - notch filt & bandpass filt EEG signals
    - convert the EEG signal from shape = (channels_num, sample_num) to (windows_num, channels_num, windows_size)
    - use 'step_time' for overlapping

    Arguments:
    EEG -- np.ndarray, size=(channel_num, sample_num)
    window_time -- the duration of each segment (window) in seconds
    step_time -- the step size (overlap) between consecutive windows in seconds. For this example, step_time==window_time, no overlap.
    Fz -- The sampling frequency (Hz) of the EEG signal

    Keyword arguments:
    notch_freq
    bandpass_freq -- a tuple (low_freq, high_freq)
    start_end_remove_window_num -- default 0, number of windows removed at the beginning and the end of the EEG signal
    amplitude_thres -- default 1000, mark all segments with np.any(EEG_seg>=amplitude_thres)=True
    to_remove_mean -- default False, whether to remove the mean of EEG signal from each channel

    Outputs:    
    EEG segments -- a list of np.ndarray of length = windows_num, each has size=(channel_num, window_size)
    segment start ids -- a list of starting ids of each window in (sample_num,)
    """
    # std_thres1 = 0.2
    # std_thres2 = 0.5
    # flat_seconds = 2
    #EEG = EEG_[[0,1,3,4]]
    #referential montage...
    
    if to_remove_mean:
        EEG = EEG - np.mean(EEG,axis=1, keepdims=True)
    window_size = int(round(window_time*Fs)) # Number of samples in each window (window_time * Fs).
    step_size = int(round(step_time*Fs)) # Step size between consecutive windows (step_time * Fs).
    # flat_length = int(round(flat_seconds*Fs))
    ## start_ids
    
    start_ids = np.arange(0, EEG.shape[1]-window_size+1, step_size) # An array of indices where each segment starts.
    if start_end_remove_window_num>0:
        # Removes start_end_remove_window_num windows from the beginning and end.
        start_ids = start_ids[start_end_remove_window_num:-start_end_remove_window_num]
    if len(start_ids) <= 0:
        raise ValueError('No EEG segments')
    
    
    ## apply montage (interchangeable with linear filtering)
        
    #EEG_segs = EEG_segs[:,[0,1,0,2],:]-EEG_segs[:,[2,3,1,3],:]
    #EEG = EEG[[0,1,0,2]] - EEG[[2,3,1,3]]
    
    ## filter signal
    #import pdb;pdb.set_trace()
    if np.max(bandpass_freq)>=notch_freq:
        EEG = notch_filter(EEG, Fs, notch_freq, n_jobs=n_jobs, verbose='ERROR')  # (#window, #ch, window_size+2padding)
    EEG = filter_data(EEG, Fs, bandpass_freq[0], bandpass_freq[1], n_jobs=n_jobs, verbose='ERROR')  # take the value starting from *padding*, (#window, #ch, window_size+2padding)
    
    EEG_segs = EEG[:,list(map(lambda x:np.arange(x,x+window_size), start_ids))].transpose(1,0,2)  # (#window, #ch, window_size+2padding)

    return EEG_segs


