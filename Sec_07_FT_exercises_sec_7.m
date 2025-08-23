%% Exercises from section 7

clear, clc

%% Exercise 1

% This question is similar to exercise #3 from the FFT section, but here using real data. The online
% material includes a file called LFPdata.mat, which contains 200 repetitions of brain electrical
% activity (LFP=local field potential), where time=0 is when a picture appeared on a computer
% screen. Compute and show the following power spectra:
%   1. Average the time-domain repetitions together, and then compute the power spectrum.
%   2. Compute the FFT for each repetition, then average the complex Fourier coefficients together,
%      and then extract the power spectrum.
%   3. Compute the FFT for each repetition, then extract the power spectrum separately for each
%      repetition, and then average the power spectra.

load LFPdata.mat

data = lfp;
pnts = length(timevec);

% (1)
spect1 = (2*abs(fft( mean(data,2) )/pnts)).^2;

% (2)
spect2 = (2*abs( mean(fft(data,[],1)/pnts ,2))).^2;

% (3)
spect3 = mean((2*abs(fft(data,[],1)/pnts)).^2, 2);


hz = linspace(0,srate/2,floor(pnts/2)+1);

% plotting
figure(1), clf
plot(hz,spect1(1:length(hz)),'s-','linew',2,'markersize',7)
hold on
plot(hz,spect2(1:length(hz)),'o-','linew',1.5,'markersize',5)
plot(hz,spect3(1:length(hz)),'p-','linew',1.5,'markersize',5)
legend({'time average';'coefs average';'power average'})
set(gca,'xlim',[0 100])
xlabel('Frequency (Hz)'), ylabel('Power')
zoom on

%% Exercise 2

% Create a multipolar chirp using frequencies that change over time in a triangular fashion. The
% following steps will help you:
%   1. Create a triangle time series. There are several ways to compute a triangle time series; one
%      equation you can use is a|mt mod 2 -1|, where a is an amplitude modulating factor, t
%      mod 2 means the modulus, or remainder, of dividing time vector t by 2, m is a frequency
%      modulator, and | | indicates the absolute value.
%      This triangle time series is the vector of frequencies at each time point.
%   2. Compute a vector y as the mean-centered cumulative sum of the triangle computed above.
%   3. Finally, the signal can be defined as s = sin(2pift + y2pi/q)
%      where f is the average of the triangle time series, t is time, and q is the sampling rate in Hz
%      (assuming t is in seconds).

% simulation details
srate = 1000;
t = 0:1/srate:5;
n = length(t);
a = 10; % amplitude modulator
m = 2; % frequency modulator
freqTS = a*abs(mod(m*t,2)-1);

% create signal
cf = mean(freqTS);
k = 2*pi/srate;
y = sin(2*pi.*cf.*t + k*cumsum(freqTS-cf));

% plotting
figure(2), clf
subplot(211), 
plot(t,freqTS,'linew',1.5)
xlabel('Time (s)'), ylabel('Frequency (Hz)')

subplot(212), 
plot(t,y,'linew',1.5)
xlabel('Time (s)'), ylabel('Amplitude')

%% Exercise 3

% Show the amplitude spectrum of the chirp. Is it interpretable? Does 
% that depend on the parameters you select (a and m are the primary parameters 
% to focus on changing)?

hz = linspace(0,srate/2,floor(n/2)+1);

amp = abs(fft(y)/n);
amp(2:end) = 2*amp(2:end);

figure(3), clf
subplot(211), 
plot(t,y,'linew',1.5)
xlabel('Time (s)'), ylabel('Amplitude')

subplot(212)
stem(hz,amp(1:length(hz)),'ks-','linew',2)
set(gca,'xlim',[0 30])
xlabel('Frequency (Hz)'), ylabel('Amplitude')

% you can manipulate the (non-)stationarity by manipulating a and m 

%% Exercise 4

% Clearly, this signal is non-stationary. Is it possible to get a more visually interpretable amplitude
% spectrum by taking the FFT of a subset of the time series?

% try Welch's method

% parameters
winlen = 1250; % window length in points (not ms!)
nbins = floor(length(t)/winlen);

% vector of frequencies for the small windows
hzL = linspace(0,srate/2,floor(winlen/2)+1);

% initialize time-frequency matrix
welchspect = zeros(1,length(hzL));

% Hann taper
hwin = .5*(1-cos(2*pi*(1:winlen) / (winlen-1)));

% loop over time windows
for ti=1:nbins
    
    % extract part of the signal
    tidx    = (ti-1)*winlen+1:ti*winlen;
    tmpdata = y(tidx);
    
    % FFT of these data (does the taper help?)
    x = fft(hwin.*tmpdata)/winlen;
    
    % and put in matrix
    welchspect = welchspect + 2*abs(x(1:length(hzL)));
end

% divide by nbins to complete average
welchspect = welchspect/nbins;

figure(4), clf
subplot(311), 
plot(t,y,'linew',1.5)
xlabel('Time (s)'), ylabel('Amplitude')
title('Signal')

subplot(312)
stem(hz,amp(1:length(hz)),'ks-','linew',2)
set(gca,'xlim',[0 30])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Frequency domain (FFT)')

subplot(313)
stem(hzL,welchspect,'ks-','linew',2)
set(gca,'xlim',[0 30],'ylim',[0 0.6])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Frequency domain (Welch''s method)')


%% Exercise 5

% I didn’t talk in much detail about the mechanisms of time-frequency analyses, but the short-time
% FFT is just a minor extension of the ”full-time” FFT. Copy/paste the code from the ”Solutions...”
% video to implement a time-frequency analysis on this signal. Is this result more interpretable? Try
% it again with different parameters for the time series.

% ST-FFT

winlen   = 500; % window length
stepsize = 25;  % step size for STFFT
numsteps = floor( (n-winlen)/stepsize );

hz = linspace(0,srate/2,floor(winlen/2)+1);


% initialize time-frequency matrix
tf = zeros(length(hz),numsteps);

% Hann taper
hwin = .5*(1-cos(2*pi*(1:winlen) / (winlen-1)));

% loop over time windows
for ti=1:numsteps
    
    % extract part of the signal
    tidx    = (ti-1)*stepsize+1:(ti-1)*stepsize+winlen;
    tapdata = y(tidx);
    
    % FFT of these data
    x = fft(hwin.*tapdata)/winlen;
    
    % and put in matrix
    tf(:,ti) = 2*abs(x(1:length(hz)));
end

figure(5), clf
subplot(311), 
plot(t,y,'linew',1.5)
xlabel('Time (s)'), ylabel('Amplitude')
title('Signal')

subplot(212)
contourf(t( (0:numsteps-1)*stepsize+1 ),hz,tf,40,'linecolor','none')
set(gca,'ylim',[0 50],'xlim',[0 5],'clim',[0 .5])
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Time-frequency power via short-time FFT')
subtitle('100 ms window')
colorbar

%% Exercise 6

% f you copy/paste the STFFT code with minimal adjustments, does the time-frequency power
% plane match the temporal dynamics of the signal? The correct answer is ”yeah sort of, but not
% exactly.” Your goal is to explain why this happens.

% there is a shift due to the witdth of the time-windows, which here
% consists of 500 ms and thus generate quite a massive averaging (try
% shorter windows)



