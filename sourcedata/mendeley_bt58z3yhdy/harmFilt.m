function [Hd,b,a]=harmFilt(r,F0,Fs,showFig)
%this is an adaptation of the Matlab comb filter, for use when the 
%sampling rate is not a multiple of the line noise frequency

w=F0*2*pi/Fs; %angle
P=floor(Fs/2/F0);
C=2*cos([1:P].*w);
B=[ones(P,1) -C' ones(P,1)];
A=[ones(P,1) -r.*C' (r^2).*ones(P,1)];
sos=[B A];
[b,a]=sos2tf(sos);
Hd = dfilt.df1tsos(sos);
if showFig
   fvtool(sos,'Fs',Fs);
   fvtool(b,a,'Fs',Fs);
end

