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

"""
This script processes individual EEG files(edf) using segment_EEG.py and segment_EEG_without_detection.py.
1. Read file and resample, filter
2. Convert EEG to Bipolar Montage
3. Segment monopolar EEG (without artifact reduction) and bipolar EEG (with artifact reduction)
4. Save processed EEG data into a MATLAB .mat file.
"""

import pdb
import datetime
import time
from collections import Counter#, deque
import os
# import os.path
import pickle
import sys
import subprocess
import scipy
#import matlab.engine
import h5py
import hdf5storage as hs
import scipy.io as sio
import numpy as np
import pandas as pd
from segment_EEG import *
from segment_EEG_without_detection import *
#from mne import *
import mne
import math


## Things to do 1/2
# If 10-20 system (MGH)
available_channels = ['EEG Fp1', 'EEG F3', 'EEG C3', 'EEG P3', 'EEG F7', 'EEG T3', 'EEG T5', 'EEG O1', 'EEG Fz', 'EEG Cz',
                      'EEG Pz', 'EEG Fp2',  'EEG F4', 'EEG C4', 'EEG P4', 'EEG F8', 'EEG T4', 'EEG T6', 'EEG O2'] # MGH BWH ULB
bipolar_channels = ['Fp1-F7', 'F7-T3', 'T3-T5', 'T5-O1', 'Fp2-F8', 'F8-T4', 'T4-T6', 'T6-O2', 'Fp1-F3', 'F3-C3', 'C3-P3',
                    'P3-O1',  'Fp2-F4', 'F4-C4', 'C4-P4', 'P4-O2', 'Fz-Cz', 'Cz-Pz']
# If 10-10 system (YAL) match 10-20
# available_channels = ['EEG Fp1', 'EEG F3', 'EEG C3', 'EEG P3', 'EEG F7', 'EEG T7', 'EEG P7', 'EEG O1', 'EEG Fz', 'EEG Cz',
#                       'EEG Pz', 'EEG Fp2',  'EEG F4', 'EEG C4', 'EEG P4', 'EEG F8', 'EEG T8', 'EEG P8', 'EEG O2'] # MGH BWH ULB
# bipolar_channels = ['Fp1-F7', 'F7-T7', 'T7-P7', 'P7-O1', 'Fp2-F8', 'F8-T8', 'T8-P8', 'P8-O2', 'Fp1-F3', 'F3-C3', 'C3-P3',
#                     'P3-O1',  'Fp2-F4', 'F4-C4', 'C4-P4', 'P4-O2', 'Fz-Cz', 'Cz-Pz']


Fs = 200.
window_time = 5  # [s]
window_step = 5  # [s]
start_end_remove_window_num = 0
amplitude_thres = 500 #500  # [uV]
line_freq = 60.  # [Hz]
bandpass_freq = [0.5, 30.]  # [Hz]
tostudy_freq = [0.5, 30.]  # [Hz]
random_state = 1

seg_mask_explanation = np.array([
    'normal',
    'NaN in EEG', #_[1,3] (append channel ids)
    'overly high/low amplitude',
    'flat signal',
    'NaN in feature',
    'NaN in spectrum',
    'overly high/low total power',
    'muscle artifact',
    'multiple assessment scores',
    'spurious spectrum',
    'fast rising decreasing',
    '1Hz artifact',])

if __name__=='__main__':

    # Things To Do 2/2:
    # Change path for input and output files
    file_path = "E:\\Zhongwei\\SAH Code Publish\\RawInput\\"
    save_path = "E:\\Zhongwei\\SAH Code Publish\\Preprocessed\\"

    #file_list = sorted(file_list, reverse=True)
    #import pdb;pdb.set_trace()
    file_list = os.listdir(file_path)
    # file_list = file_list[727:]
    
    for ifile in file_list:
        file = file_path + ifile # + "\\" + ifile + "_ARed.edf"
        print(file)
        #import pdb;pdb.set_trace()
#        if os.path.isfile(save_path+ifile+'.mat'):
#            continue
#        else:
#            try:
            #data = mne.io.read_raw_edf(file,preload=True)
        data = mne.io.read_raw_edf(file,stim_channel=None,exclude='EDF Annotations',preload=True)
        raw_data = data.get_data(picks=range(23))
        #import pdb;pdb.set_trace()
        info = data.info
        fs = info['sfreq']
        #raw_data = scipy.signal.resample(raw_data, int(math.floor(raw_data.shape[1]*Fs/fs)),axis=1)
        if fs!=Fs:
            raw_data = scipy.signal.resample_poly(raw_data, Fs, fs, axis=1)
            #raw_data = mne.filter.resample(raw_data, down=fs/Fs, npad='auto')
        raw_data = raw_data*10e5 # V->uV
        
        channels = data.ch_names
        channels = [x.upper() for x in channels]
        chan_index = list()
        for chNo in available_channels:
            chan_index.append(channels.index(chNo.upper()))
            # chan_index.append(channels.index(chNo.upper() + '-REF'))
        raw_data = raw_data[chan_index,:]
        
        
        
        ## Bipolar reference
        # The EEG signal is computed as the voltage difference between two adjacent electrodes,
        # which helps reduce common noise artifacts
        # monopolar: referenced to a single common reference
        bipolar_data = np.zeros((18,raw_data.shape[1]))
        bipolar_data[8,:] = raw_data[0,:] - raw_data[1,:]; # Fp1-F3
        bipolar_data[9,:] = raw_data[1,:] - raw_data[2,:]; # F3-C3
        bipolar_data[10,:] = raw_data[2,:] - raw_data[3,:]; # C3-P3
        bipolar_data[11,:] = raw_data[3,:] - raw_data[7,:]; # P3-O1
    
        bipolar_data[12,:] = raw_data[11,:] - raw_data[12,:]; # Fp2-F4
        bipolar_data[13,:] = raw_data[12,:] - raw_data[13,:]; # F4-C4
        bipolar_data[14,:] = raw_data[13,:] - raw_data[14,:]; # C4-P4
        bipolar_data[15,:] = raw_data[14,:] - raw_data[18,:]; # P4-O2
    
        bipolar_data[0,:] = raw_data[0,:] - raw_data[4,:];  # Fp1-F7
        bipolar_data[1,:] = raw_data[4,:] - raw_data[5,:]; # F7-T3
        bipolar_data[2,:] = raw_data[5,:] - raw_data[6,:]; # T3-T5
        bipolar_data[3,:] = raw_data[6,:] - raw_data[7,:]; # T5-O1
    
        bipolar_data[4,:] = raw_data[11,:] - raw_data[15,:]; # Fp2-F8
        bipolar_data[5,:] = raw_data[15,:] - raw_data[16,:]; # F8-T4
        bipolar_data[6,:] = raw_data[16,:] - raw_data[17,:]; # T4-T6
        bipolar_data[7,:] = raw_data[17,:] - raw_data[18,:]; # T6-O2
    
        bipolar_data[16,:] = raw_data[8,:] - raw_data[9,:];   # Fz-Cz
        bipolar_data[17,:] = raw_data[9,:] - raw_data[10,:]; # Cz-Pz
        
        ## save 5s monopolar/bipolar epoches using notch/band pass/artifact detection/resampling

        segs_monpolar = segment_EEG_without_detection(raw_data,available_channels,window_time, window_step, Fs,
                            notch_freq=line_freq, bandpass_freq=bandpass_freq,
                            to_remove_mean=False, amplitude_thres=amplitude_thres, n_jobs=-1, start_end_remove_window_num=start_end_remove_window_num)
        del raw_data
        segs_, bs_, seg_start_ids_, seg_mask, specs_, freqs_ = segment_EEG(bipolar_data,bipolar_channels,window_time, window_step, Fs,
                            notch_freq=line_freq, bandpass_freq=bandpass_freq,
                            to_remove_mean=False, amplitude_thres=amplitude_thres, n_jobs=-1, start_end_remove_window_num=start_end_remove_window_num)
        # segs_ = segment_EEG_without_detection(bipolar_data, bipolar_channels, window_time,
        #                                                                    window_step, Fs,
        #                                                                    notch_freq=line_freq,
        #                                                                    bandpass_freq=bandpass_freq,
        #                                                                    to_remove_mean=False,
        #                                                                    amplitude_thres=amplitude_thres, n_jobs=-1,
        #                                                                    start_end_remove_window_num=start_end_remove_window_num)

        if len(segs_) <= 0:
            raise ValueError('No segments')
            
        seg_mask2 = map(lambda x:x.split('_')[0], seg_mask)
        sm = Counter(seg_mask2)
        for ex in seg_mask_explanation:
            if ex in sm:
                print('%s: %d/%d, %g%%'%(ex,sm[ex],len(seg_mask),sm[ex]*100./len(seg_mask)))
                
        if segs_.shape[0]<=0:
            raise ValueError('No EEG signal')
            if segs_.shape[1]!=len(bipolar_channels):
                raise ValueError('Incorrect #chanels')
        
        fd = os.path.split(save_path)[0]
        if not os.path.exists(fd):
            os.mkdir(fd)
        res = {'EEG_segs_bipolar':segs_.astype('float16'),
               'EEG_segs_monopolar':segs_monpolar.astype('float16'),
               'EEG_specs':specs_.astype('float16'),
               'burst_suppression':bs_.astype('float16'),
               'EEG_frequency':freqs_,
               'seg_start_ids':seg_start_ids_,
               'Fs':Fs,
               'seg_masks':seg_mask,
               'channel_names':bipolar_channels}
        # res = {'EEG_segs_bipolar': segs_.astype('float16'),
        #        'EEG_segs_monopolar': segs_monpolar.astype('float16'),
        #        'channel_names': bipolar_channels}
        # sio.savemat(save_path+ifile +'_raw.mat', res, do_compression=True)
        sio.savemat(save_path + ifile + '.mat', res, do_compression=True)
#