function y = single_waveform_feature_demo(x, bakI, fs)
%% compute wavelet entropy feature of a single waveform with fixed length
%% x: a short(we use 200ms) intracranial EEG waveform centered around an single HFO
%% detection; 
%% bakI: subscripts of the signal of interest( we use a 100ms-long waveform center around
%% an HFO detection) 
%% fs: sampling frequency ( we use 4096 Hz)
    fsThre = [ 80 500 ]; % frequency band of interest
    if size(x,1) > 1
        x = x';
    end
    x = x-mean(x);
    [cMat,f]=cwt(x,fs); % get a time-frequency map of the signal using continuous wavelet transform
 
    % use only frequency components that are within the frequency band specified
    % in fsThre; if the upper limit in fsThre exceeds the maximum
    % frequency component contained in the signal, the later is
    % used as the upper limit.
    If = f >= fsThre(1) & f <= 0.4*fs & f <= fsThre(2);
                                                       
    a=abs(cMat(If,bakI)); % extract entries within the specified
                          % frequency/time range
    a=reshape(a,1,[]);
    a = abs(a)./sum(abs(a)); % construct a discrete probability
                             % distribution for the entries in the time-frequency map
    wEntropy = -sum(a.*(log(eps+a))); % compute Shannon entropy of the
                                      % probability distribution  
    y = wEntropy; % wavelet entropy is used as the feature to
                  % distinguish EMG-false-HFOs from real ones 

end
