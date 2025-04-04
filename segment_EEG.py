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
"""


from collections import Counter
import numpy as np
from joblib import Parallel, delayed
from scipy.signal import detrend
import matplotlib.pyplot as plt
import mne
mne.set_log_level(verbose='WARNING')
from mne.filter import filter_data, notch_filter
from peakdetect import peakdetect
#from read_delirium_data import datenum


seg_mask_explanation = [
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
    '1Hz artifact',]
    


# Artifact Detection Function:

def peak_detect_num_amp(signal, lookahead=200, delta=0):

    # Artifact Detection Functions:
    # Detect peaks in EEG siganls and returns the number of peaks and their amplitude range.
    # signal: #channel x #points
    res_num = []
    res_amp = []
    for cid in range(signal.shape[0]):
        local_max, local_min = peakdetect(signal[cid], lookahead=lookahead, delta=delta)
        if len(local_min)<=0:
            local_extremes = np.array(local_max)
        elif len(local_max)<=0:
            local_extremes = np.array(local_min)
        else:
            local_extremes = np.r_[local_max, local_min]
        res_num.append(len(local_extremes))
        if len(local_extremes)<=0:
            amp = 0
        else:
            amp = local_extremes[:,1].max()-local_extremes[:,1].min()
        res_amp.append(amp)
    return res_num, res_amp
    
    
def peak_detect(signal, max_change_points, min_change_amp, lookahead=200, delta=0):

    # Artifact Detection Function:
    # Detects fast changes in amplitude (e.g., muscle artifacts).
    # signal: #channel x #points
    res = []
    for cid in range(signal.shape[0]):
        local_max, local_min = peakdetect(signal[cid], lookahead=lookahead, delta=delta)
        if len(local_min)<=0 and len(local_max)<=0:
            res.append(False)
        else:
            if len(local_min)<=0:
                local_extremes = np.array(local_max)
            elif len(local_max)<=0:
                local_extremes = np.array(local_min)
            else:
                local_extremes = np.r_[local_max, local_min]
            local_extremes = local_extremes[np.argsort(local_extremes[:,0])]
            res.append(np.logical_and(np.diff(local_extremes[:,0])<=max_change_points, np.abs(np.diff(local_extremes[:,1]))>=min_change_amp).sum())
    return res
    
    
def smooth(x, window_len=12, window='hanning'):

    # Artifact Detection Function:
    # Applies a moving average filter to smooth the EEG signal.
    #'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
    s = np.r_[x[window_len//2-1:0:-1],x,x[-2:-window_len//2-1:-1]]
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = eval('np.'+window+'(window_len)')

    return np.convolve(w/w.sum(), s, mode='valid')


def autocorrelate_noncentral_max_abs(x):

    # Artifact Detection Function:
    # Computes the autocorrelation to detect repetitive signal patterns (useful for identifying ECG artifacts in EEG).
    ress = []
    for ii in range(x.shape[1]):
        res = np.correlate(x[:,ii],x[:,ii],mode='full')/np.correlate(x[:,ii],x[:,ii],mode='valid')[0]
        ress.append(np.max(res[len(res)//2+7:len(res)//2+20]))  # ECG range: 40/1min(0.7Hz) -- 120/1min(2Hz)
    return ress





def segment_EEG(EEG, Ch_names, window_time, step_time, Fs, notch_freq=None, bandpass_freq=None, start_end_remove_window_num=0,
                amplitude_thres=500, n_jobs=1, to_remove_mean=False):#
    """Segment EEG signals with artifact detection.
    - step1. Mean Removal (if enabled)
    - step2. Segment EEG into Windows (overlapped depending on step_time)
    - step3. Filtering
    - step4. Detect and Remove Artifacts
            - NaN checks
            - Large amplitude fluctuations
            - Flat EEG signals (low variance)
            - Fast rising/falling signals (artifacts)
            - ECG contamination detection (based on autocorrelation).
    - step5. Compute Spectrograms for EEG Analysis

    Arguments:
    EEG -- Input EEG signals, np.ndarray(channel_num, sample_num)
    Ch_names -- List of EEG channel names.
    window_time -- Length of each segment in seconds.
    step_time -- Step size for segmentation.
    Fz -- Sampling frequency (Hz).
    notch_freq -- Frequency for notch filtering (removes power line noise).
    bandpass_freq -- Bandpass filter range (low and high cutoff frequencies).
    start_end_remove_window_num -- default 0, number of windows removed at the beginning and the end of the EEG signal
    amplitude_thres -- Threshold for detecting large amplitude artifacts, default 1000, mark all segments with np.any(EEG_seg>=amplitude_thres)=True
    to_remove_mean -- default False, whether to remove the mean of EEG signal from each channel

    Outputs:    
    EEG segments -- a list of np.ndarray of length = windows_num, each has size=(channel_num, window_size)
    labels --  a list of labels of each window
    segment start ids -- a list of starting ids of each window in (sample_num,)
    segment masks --
    """
    std_thres1 = 0.2
    std_thres2 = 0.5
    flat_seconds = 2
    #EEG = EEG_[[0,1,3,4]]
    #referential montage...

    # segmentation step1 - Mean Removal (if enabled):
    if to_remove_mean:
        EEG = EEG - np.mean(EEG,axis=1, keepdims=True)

    window_size = int(round(window_time*Fs))
    step_size = int(round(step_time*Fs))
    flat_length = int(round(flat_seconds*Fs))
    ## start_ids

    # segmentation step2 - Segment EEG into Windows
    start_ids = np.arange(0, EEG.shape[1]-window_size+1, step_size)
    if start_end_remove_window_num>0:
        start_ids = start_ids[start_end_remove_window_num:-start_end_remove_window_num]
    if len(start_ids) <= 0:
        raise ValueError('No EEG segments')
    
    seg_masks = [seg_mask_explanation[0]]*len(start_ids)
    """
    for i in range(len(start_ids)):
        ll = labels[start_ids[i]:start_ids[i]+window_size]
        nll = np.isnan(ll)
        if not np.all(nll) and (np.any(nll) or len(set(ll))!=1):
            seg_masks[i] = seg_mask_explanation[8]
    """
    
    ## apply montage (interchangeable with linear filtering)
        
    #EEG_segs = EEG_segs[:,[0,1,0,2],:]-EEG_segs[:,[2,3,1,3],:]
    #EEG = EEG[[0,1,0,2]] - EEG[[2,3,1,3]]
    
    ## segmentation step3 - filter signal
    #import pdb;pdb.set_trace()
    if np.max(bandpass_freq)>=notch_freq:
        EEG = notch_filter(EEG, Fs, notch_freq, n_jobs=n_jobs, verbose='ERROR')  # (#window, #ch, window_size+2padding)
    EEG = filter_data(EEG, Fs, bandpass_freq[0], bandpass_freq[1], n_jobs=n_jobs, verbose='ERROR')  # take the value starting from *padding*, (#window, #ch, window_size+2padding)
    
    ## detect burst suppression
    # BSR_segs is related to burst suppression ratio (BSR).It quantifies the fraction of time the EEG shows burst suppression
    # (periods of high activity followed by suppression).
    # Computed from the BS (burst suppression envelope data) and helps quantify the presence of burst suppression over time windows.
    # BSR is a biomarker for brain function, especially in:
    # Coma patients
    # Anesthesia monitoring
    # Neurocritical care
    # A high BSR value indicates a longer suppression period (bad prognosis).
    # A low BSR value suggests more continuous EEG activity (good prognosis).
    

    # 1.Compute the Hilbert Envelope of the EEG Signal

    # applies smoothing to the EEG signal before processing, Uses a moving average filter with a window_len=10.
    EEG_tmp = np.zeros_like(EEG)
    for i in range(EEG.shape[0]):
        eeg_smooth = smooth(EEG[i,:], window_len=10, window='flat')
        EEG_tmp[i,:eeg_smooth.shape[0]] = eeg_smooth

    # Converts EEG data into an MNE RawArray.
    EEG_mne = mne.io.RawArray(np.array(EEG_tmp, copy=True), mne.create_info(Ch_names, Fs, ch_types='eeg', verbose='ERROR'), verbose='ERROR')
    # Applies the Hilbert Transform to compute the EEG envelope (absolute value of the analytic signal).
    EEG_mne.apply_hilbert(envelope=True, n_jobs=-1, verbose='ERROR')
    BS = EEG_mne.get_data()


    # 2.Define the Burst Suppression Window
    bs_window_size = int(round(120*Fs))  # BSR is computed over a 2min window
    bs_start_ids = np.arange(0, EEG.shape[1]-bs_window_size+1,bs_window_size) # starting indices of the 2-minute windows.

    # Handling edge cases: ensures that even if there isn’t an exact 2-minute segment, the last window still gets processed.
    if len(bs_start_ids)<=0:
        bs_start_ids = np.array([0], dtype=int)
    if EEG.shape[1]>bs_start_ids[-1]+bs_window_size:  # if incomplete divide (if there’s leftover data at the end
        bs_start_ids = np.r_[bs_start_ids, EEG.shape[1]-bs_window_size]


    # 3.Compute Burst Suppression Ratio (BSR_segs)
    BS_segs = BS[:,list(map(lambda x:np.arange(x,min(BS.shape[1],x+bs_window_size)), bs_start_ids))] # non-overlapping windows
    BSR_segs = np.sum(BS_segs<=5, axis=2).T*1./bs_window_size
    #Computes the fraction of time where the Hilbert envelope ≤ 5 µV (indicating suppression).
    #Result normalized by the segment length (bs_window_size) to get the percentage of suppression time

    # 4.Assign BSR to the Entire EEG Signal, BSR is mapped to the entire EEG time series.
    BSR = np.zeros_like(EEG)
    for ii, bsi in enumerate(bs_start_ids):
        BSR[:, bsi:min(BSR.shape[1], bsi+bs_window_size)] = BSR_segs[ii].reshape(-1,1)


    # 5.Segment BSR into EEG Windows
    # Segments BSR values into smaller EEG time windows (defined by window_size).
    # Averaging over each window gives a single BSR value per EEG window.
    BSR_segs = BSR[:,list(map(lambda x:np.arange(x,x+window_size), start_ids))].transpose(1,0,2).mean(axis=2)
    EEG_segs = EEG[:,list(map(lambda x:np.arange(x,x+window_size), start_ids))].transpose(1,0,2)  # (#window, #ch, window_size+2padding)
    #EEG_segs = detrend(EEG_segs, axis=2)
    #EEG_segs = EEG[:,map(lambda x:np.arange(x-padding,x+window_size+padding), start_ids)].transpose(1,0,2)  # (#window, #ch, window_size+2padding)
    #if np.max(bandpass_freq)>=notch_freq:
    #    EEG_segs = notch_filter(EEG_segs, Fs, notch_freq, n_jobs=n_jobs, verbose='ERROR')  # (#window, #ch, window_size+2padding)
    #EEG_segs = filter_data(detrend(EEG_segs, axis=2), Fs, bandpass_freq[0], bandpass_freq[1], n_jobs=n_jobs, verbose='ERROR')  # take the value starting from *padding*, (#window, #ch, window_size+2padding)
    #EEG_segs = EEG_segs[:,:,padding:-padding]  # (#window, #ch, window_size)detrend(, axis=2)
    
    ## find NaN in signal
    # If a segment contains NaN values in EEG, it is labeled accordingly
    nan2d = np.any(np.isnan(EEG_segs), axis=2)
    nan1d = np.where(np.any(nan2d, axis=1))[0]
    for i in nan1d:
        seg_masks[i] = '%s_%s'%(seg_mask_explanation[1], np.where(nan2d[i])[0])
        
    ## calculate spectrogram
    # Computes power spectral density (PSD): Uses multitaper method with a bandwidth of 2 Hz to estimate the power in different frequency bands.
    # Converts power to dB: Transposes and converts values into logarithmic scale (10 * log10(power)).
    # Defines frequency resolution (df): The spacing between frequency bins.
    
    #TODO detrend(EEG_segs)
    #TODO remove_mean(EEG_segs) to remove frequency at 0Hz
    
    #mne_epochs = mne.EpochsArray(EEG_segs, mne.create_info(ch_names=list(map(str, range(EEG_segs.shape[1]))), sfreq=Fs, ch_types='eeg'), verbose=False)
    BW = 2.
    specs, freq = mne.time_frequency.psd_array_multitaper(EEG_segs, Fs, fmin=bandpass_freq[0], fmax=bandpass_freq[1], adaptive=False, low_bias=False, n_jobs=n_jobs, verbose='ERROR', bandwidth=BW, normalization='full')
    df = freq[1]-freq[0]
    specs = 10*np.log10(specs.transpose(0,2,1))
    
    ## find nan in spectrum
    # Replaces inf values with NaN to handle log-transform issues.
    # Finds spectrograms with NaNs and logs them.
    
    specs[np.isinf(specs)] = np.nan
    nan2d = np.any(np.isnan(specs), axis=1)
    nan1d = np.where(np.any(nan2d, axis=1))[0]
    nonan_spec_id = np.where(np.all(np.logical_not(np.isnan(specs)), axis=(1,2)))[0]
    for i in nan1d:
        seg_masks[i] = '%s_%s'%(seg_mask_explanation[5],np.where(nan2d[i])[0])
        
    ## find staircase-like spectrum
    # | \      +-+
    # |  \     | |
    # |   -----+ +--\
    # +--------------=====
    # Detects "staircase-like" artifacts in EEG spectrograms.
    # Uses two filters (aa for increasing patterns, bb for decreasing patterns) to identify abrupt spectral changes.
    # If detected, the affected segment indices are logged.

    spec_smooth_window = int(round(1./df))  # 1 Hz
    specs2 = specs[nonan_spec_id][:,np.logical_and(freq>=5,freq<=20)]
    freq2 = freq[np.logical_and(freq>=5,freq<=20)][spec_smooth_window:-spec_smooth_window]
    ww = np.hanning(spec_smooth_window*2+1)
    ww = ww/ww.sum()
    smooth_specs = np.apply_along_axis(lambda m: np.convolve(m, ww, mode='valid'), axis=1, arr=specs2)
    dspecs = specs2[:,spec_smooth_window:-spec_smooth_window]-smooth_specs
    #dspecs_std = np.std(dspecs, axis=1, keepdims=True)
    #dspecs_std[dspecs_std<1e-3] = 1.
    dspecs = dspecs-dspecs.mean(axis=1,keepdims=True)#()/dspecs_std
    aa = np.apply_along_axis(lambda m: np.convolve(m, np.array([-1.,-1.,0,1.,1.,1.,1.]), mode='same'), axis=1, arr=dspecs)  # increasing staircase-like pattern
    bb = np.apply_along_axis(lambda m: np.convolve(m, np.array([1.,1.,1.,1.,0.,-1.,-1.]), mode='same'), axis=1, arr=dspecs)  # decreasing staircase-like pattern
    stsp2d = np.logical_or(np.maximum(aa,bb).max(axis=1)>=10, np.any(np.abs(np.diff(specs2,axis=1))>=11, axis=1))
    #stsp2d = np.logical_and(np.maximum(aa,bb)>=14, np.abs(np.concatenate([dspecs[:,2:]-dspecs[:,:-2], np.zeros((dspecs.shape[0],2,dspecs.shape[2]))))>=1.8, axis=1))
    #stsp2d = np.any(np.logical_and(np.concatenate([np.zeros((dspecs.shape[0],1,dspecs.shape[2])),np.abs(np.diff(np.arctan(np.diff(dspecs,axis=1)/df),axis=1))/np.pi*180],axis=1)[:,:-1]>80, np.logical_or(np.abs(np.diff(dspecs,axis=1))[:,:-1]>2,np.abs(dspecs[:,2:]-dspecs[:,:-2])>3)),axis=1)
    #stsp2d = np.logical_or(stsp2d, np.any(np.abs(np.diff(specs2,axis=1))>=6, axis=1))
    stsp1d = nonan_spec_id[np.any(stsp2d, axis=1)]
    """
    stsp2d = Parallel(n_jobs=n_jobs,verbose=True)(delayed(peak_detect_num_amp)(dspecs[sid].T, lookahead=5, delta=1.5) for sid in range(dspecs.shape[0]))
    stsp2d = (np.array(map(lambda x:x[0],stsp2d))>=8) & (np.array(map(lambda x:x[1],stsp2d))>=6)
    stsp1d = nonan_spec_id[np.where(np.any(stsp2d, axis=1))[0]]
    """
    for i in stsp1d:
        seg_masks[i] = seg_mask_explanation[9]
    
    ## check ECG in spectrum (~1Hz and harmonics)
    #dspecs2 = dspecs[:,np.logical_and(freq2>=6,freq2<=10)]
    autocorrelation = Parallel(n_jobs=n_jobs,prefer="threads",verbose=True)(delayed(autocorrelate_noncentral_max_abs)(spec) for spec in dspecs)
    autocorrelation = np.array(autocorrelation)
    ecg2d = autocorrelation>0.7
    ecg1d = nonan_spec_id[np.any(ecg2d,axis=1)]
    for i in ecg1d:
        seg_masks[i] = seg_mask_explanation[11]


    ## find overly fast rising/decreasing signal
    # Detects sudden increases/decreases in EEG amplitude using peak_detect().
    max_change_points = 0.1*Fs
    min_change_amp = 1.8*amplitude_thres
    fast_rising2d = Parallel(n_jobs=n_jobs,prefer="threads",verbose=True)(delayed(peak_detect)(EEG_segs[sid], max_change_points, min_change_amp, lookahead=50, delta=0) for sid in range(EEG_segs.shape[0]))
    fast_rising2d = np.array(fast_rising2d)>0
    fast_rising1d = np.where(np.any(fast_rising2d, axis=1))[0]
    for i in fast_rising1d:
        seg_masks[i] = seg_mask_explanation[10]


    ## find large amplitude in signal
    # Finds extreme amplitude variations in EEG signals and logs them.
    amplitude_large2d = np.max(EEG_segs,axis=2)-np.min(EEG_segs,axis=2)>2*amplitude_thres
    amplitude_large1d = np.where(np.any(amplitude_large2d, axis=1))[0]
    for i in amplitude_large1d:
        seg_masks[i] = '%s_%s'%(seg_mask_explanation[2], np.where(amplitude_large2d[i])[0])


    ## find flat signal
    # Finds segments where EEG remains nearly constant (flat-line detection).
    # Detects brain suppression states (e.g., burst suppression).
    # careful about burst suppression
    EEG_segs_temp = EEG_segs[:,:,:(EEG_segs.shape[2]//flat_length)*flat_length]
    short_segs = EEG_segs_temp.reshape(EEG_segs_temp.shape[0], EEG_segs_temp.shape[1], EEG_segs_temp.shape[2]//flat_length, flat_length)
    flat2d = np.any(detrend(short_segs, axis=3).std(axis=3)<=std_thres1, axis=2)
    flat2d = np.logical_or(flat2d, np.std(EEG_segs,axis=2)<=std_thres2)
    flat1d = np.where(np.any(flat2d, axis=1))[0]
    for i in flat1d:
        seg_masks[i] = '%s_%s'%(seg_mask_explanation[3], np.where(flat2d[i])[0])
    
    """
    ## local z-transform
    eeg = np.array(EEG_segs[:,:,:100], copy=True)
    eeg[np.isinf(eeg)] = np.nan
    mu = np.nanmean(eeg, axis=2)
    sigma2 = np.nanstd(eeg, axis=2)**2
    EEG_segs_ = np.array(EEG_segs, copy=True)
    for t in range(EEG_segs.shape[2]):
        eeg = np.array(EEG_segs[:,:,t], copy=True)
        eeg[np.logical_or(np.isnan(eeg),np.isinf(eeg))] = 0.
        
        sigma2_no0 = np.array(sigma2, copy=True)
        sigma2_no0[np.abs(sigma2_no0)<=1e-3] = 1.
        
        eeg = (eeg-mu)/np.sqrt(sigma2_no0)
        mu = 0.001*eeg+0.999*mu
        sigma2 = 0.001*(eeg-mu)**2+0.999*sigma2
        
        EEG_segs_[:,:,t] = eeg
    EEG_segs = detrend(EEG_segs_, axis=2)
    """

    ## Computes final spectrogram with finer frequency resolution (1 Hz).

    BW = 1.  #frequency resolution 1Hz

    # Computes spectral power across frequency bands.
    specs, freq = mne.time_frequency.psd_array_multitaper(EEG_segs, Fs, fmin=bandpass_freq[0], fmax=bandpass_freq[1], adaptive=False, low_bias=False, n_jobs=n_jobs, verbose='ERROR', bandwidth=BW, normalization='full')
    specs = 10*np.log10(specs.transpose(0,2,1))
    
            
    return EEG_segs, BSR_segs, start_ids, seg_masks, specs, freq


