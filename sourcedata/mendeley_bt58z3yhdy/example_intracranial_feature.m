clear all
iEEGfileName = './intracranialEEG.mat'; 
hfoFileName = './qHFOdetection.mat'; 
 
fs = 5000; %% sampling frequency of raw iEEG
deciFs = 4096; %% sampling frequency used when computing feature
filtRange = [80 500]; %% frequency band of interest
secWin = 0.2; % length(in second) of signals used to construct
              % time-frequency map around individual HFO detections via wavelet transformation
secBak = 0.1;  % length(in second) of signals used compute
                 % wavelet entropy around individual HFO detections
secBUFFER = 1; % length(in second) of the buffer used to
                      % mitigate ringing artifacts near the edges
                      % of the windowed signals

%% compute wavelet entropy of all HFOs in file hfoFileName
features = compute_feature_demo( iEEGfileName, hfoFileName, fs, deciFs,filtRange,secWin,secBak,secBUFFER);

Threshold = 8.956; %% threshold for distinguish EMG from HFO

isEMG = features>Threshold; %% large wavelet entropy indicate EMG-false-HFO

%% visualizing results
N = 1; %%  Change this to the channel of interest
hfo = load(hfoFileName);
ieeg = matfile(iEEGfileName);
ieeg = double(ieeg.data(N,:))-mean(ieeg.data);
I = hfo.chanIdx==N;
if sum(I)==0
  disp('there is no qHFO detection in this channel');
else
  startIdx = hfo.startIdx(I);
  stopIdx = hfo.stopIdx(I);
  isEMG = isEMG(I);
  isEMGlabel = {true,false};
  isEMGcolor = {'k','r'};   %will plot true qHFO in red, EMG-false-HFO in black
  plot(ieeg)
  hold on
  for i=1:2
    plot(double(startIdx(isEMG==isEMGlabel{i})),max(ieeg),[isEMGcolor{i} 'o'])
    plot(double(stopIdx(isEMG==isEMGlabel{i})),max(ieeg),[isEMGcolor{i} '+'])
  end
  hold off
end
