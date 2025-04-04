
function [out, header] = Feature_Subset( x, i, Fs )
%% IN THIS CASE X is a matrix of the EEGs

warning('off','all')
%% FEATURES
header = {'Shannon'
    'AlphaDeltaRatio'
    'TotalPow'
    'DeltaPow'
    'ThetaPow'
    'AlphaPow'
    'PerAlphaVar'
    };
featurenum = length(header);
out = zeros(1,featurenum);
%% COMPUTE THE SHANNON ENTROPY

    bin_min = -200; bin_max = 200; binWidth = 2;
    [out(1), prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    
        %% ALPHA TO DELTA RATIO.  - 0.35
     
    window=Fs*2;
    noverlap=32;
    h = spectrum.welch('Hamming',window);
    
    try
    hpsd = psd(h,x(:,i),'Fs',Fs);
    Pw = hpsd.Data;
    Fw = hpsd.Frequencies;
    out(2) = bandpower(x(:,i),Fs,[8,13])/bandpower(x(:,i),Fs,[0.5,4]);
    catch
        out(2) = 0;
    end
    
    
        %% BAND POWERS - 
    out(3) = bandpower(x(:,i),Fs,[0.5,15]);
    out(4) = bandpower(x(:,i),Fs,[0.5,4]);
    out(5) = bandpower(x(:,i),Fs,[4,7]);
    out(6) = bandpower(x(:,i),Fs,[8,15]);
    out(7) = out(6)/(out(4)+out(5)+out(6)); % The proportion of Alpha power relative to the sum of Delta, Theta, and Alpha.
end
