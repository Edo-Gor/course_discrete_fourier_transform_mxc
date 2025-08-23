%% Exercises from section 6

clear, clc

%% Exercise 1

% Create a 2-second sine wave at 11 Hz. Compute and plot its amplitude spectrum. Next, recompute
% the Fourier transform after zero-padding to a total of 4 seconds. Plot the amplitude spectrum
% after (1) dividing the Fourier coefficients by the number of time points in the original sine wave
% and (2) dividing the coefficients by the number of time points including the padded zeros. Which
% normalization returns the accurate amplitude?

srate = 1000;
time  = 0:1/srate:2;
pnts  = length(time);

signal = 2*sin(2*pi*11*time);

Fcoeffs = fft(signal,2*pnts);

Normalised1 = 2*abs(Fcoeffs) / (2*pnts);
Normalised2 = 2*abs(Fcoeffs) / pnts;

hz = linspace(0,srate/2,floor(2*pnts/2)+1);

% plotting 
figure(1), clf
subplot(211)
stem(hz,Normalised1(1:length(hz)),'linew',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Frequency spectrum (wrong normalisation)')
subtitle('2 s 11 Hz sinewave zero-padded to 4s')
set(gca,'xlim',[0 30])

subplot(212)
stem(hz,Normalised2(1:length(hz)),'linew',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Frequency spectrum (correct normalisation)')
subtitle('2 s 11 Hz sinewave zero-padded to 4s')
set(gca,'xlim',[0 30])

% only the second returns the correct amplitude because the normalisation
% must be applied with the length of the original signal, not the zero-padded
% signal length

%% Exercise 2

% The first, because the reconstruted fequencies are discrete, not
% continuous

%% Exercise 3

% The second, because it's easier to understand

%% Exercise 4

% No you cannot, the srate define the Nyquist for instance, how can you
% define it if the srate is variable ?



