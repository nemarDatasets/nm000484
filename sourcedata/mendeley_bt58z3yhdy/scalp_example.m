clear;

%% get some data
load('data.mat')
FS = 5000;

%%load parameter for the detector, there're muliple choices (cases) with
%%different sensitivity/specifity
par = load_par(2);

%%construct detector class
obj = EMGdetClass( par );

%%parameters related to seperating data into epochs
baseline = [];
epochSec = 600;
buffSize = 0;

%% run detector
obj.doAll(data,FS,baseline,epochSec,buffSize);

%% get information of EMG detections
startIdx = obj.startIdx;
stopIdx = obj.stopIdx;
chanIdx = obj.chanIdx;


%% view results for channel 'N'

N=1;  %to view channel 1
data = obj.X(N,:);   % obj.X is the same scalp EEG after common average referencing
plot(data)
hold on
plot(obj.startIdxDe(chanIdx==N),max(data(N,:)),'o')
plot(obj.stopIdxDe(chanIdx==N),max(data(N,:)),'+')
hold off
