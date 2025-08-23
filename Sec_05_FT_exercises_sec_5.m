%% Exercises from section 5

clear, clc

%% Exercise 1

% One method of scrambling a signal is to shuffle the phases of the Fourier coefficients and then
% reconstruct the time-domain data. In this exercise, generate a linear chirp and compute its inverse
% after shuffling the phases (do not alter the power values). Plot the time domain signals and their
% spectra. Is the scrambled time series recognizable as having been a chirp?

srate= 1000;
t = 0:1/srate:6;
N = length(t);
f = [1 5];
ff = linspace(f(1),mean(f),N);
data = sin(2*pi.*ff.*t);

dataX = fft(data);

phases = angle([dataX(10:end) dataX(1:9)]);
shuffdata = abs(dataX).*exp(1i*phases);

newdata = ifft(shuffdata);

% plotting 
figure(1), clf
subplot(211)
plot(t,data),
hold on
plot(t,real(newdata),'r')
legend({'Original signal','Phase-shuffled signal'})
title('Time domain')

subplot(212)
hz = linspace(0,srate/2,floor(N/2)+1);
plot(hz,2*abs(dataX(1:length(hz))/N),'o'),
hold on
plot(hz,2*abs(shuffdata(1:length(hz))/N),'r')
legend({'Original signal','Phase-shuffled signal'})
title('Frequency domain')
set(gca,'xlim',[0 20])

%% Exercise 2

% Real data often contain noise. One way to reduce noise is to average over more data; if the noise
% is random, then it will decrease with increased averaging. This and the next several exercises will
% explore this idea.
%   > Generate a signal comprising three sine waves at three frequencies with three different amplitudes. 
%   You can pick the parameters, but make sure the signal contains only peaks at the
%   frequencies, and not at neighboring frequencies due to a mismatch in the frequencies, time
%   ranges, and frequency resolution.
%   > Generate 50 noisy repetitions of the signal by creating a new instance of the signal and adding
%   random noise. It’s most convenient if you store the data in a 50 by time matrix. Make sure
%   to use a sufficient amplitude of noise so the signals look noisy.

srate = 1/1000;
time = 0:srate:2;
pnts = length(time);

freq1 = 5; freq2 = 12; freq3 = 19;
ampl1 = 1; ampl2 = 3;  ampl3 = 5;

signal = ampl1*sin(2*pi*freq1*time) + ...
         ampl2*sin(2*pi*freq2*time) + ...
         ampl3*sin(2*pi*freq3*time);

data = zeros(50,pnts);
for i = 1:50
    
    data(i,:) = signal + 50*randn(size(signal));
    
end

% plotting
figure(2), clf
subplot(211)
plot(time,signal)
title('Original signal')

subplot(212)
plot(time,data)
title('Noisy data')

%% Exercise 3

% Show three power spectra from these data: (1) average over all repetitions in the time domain
% and then computing one FFT; (2) compute the FFT of each repetition (thus 50 FFTs, note that
% you can do this by inputting a matrix into the fft function, just make sure the FFT is computed
% along the time dimension!), average the coefficients together and then extract power; (3) same
% as #2 except first extract power for each repetition and then average the power spectra together.
% Comment on any differences in the results.

% (1)
spect1 = 2*abs(fft( mean(data,1) )/pnts).^2;

% (2)
spect2 = 2*abs( mean(fft(data,[],2)/pnts ,1)).^2;

% (3)
spect3 = mean( 2*abs( fft(data,[],2)/pnts ).^2,1);

hz = linspace(0,(1/srate)/2,floor(pnts/2)+1);

% plotting
figure(3), clf
plot(hz,spect1(1:length(hz)),'s-','linew',2,'markersize',10)
hold on
plot(hz,spect2(1:length(hz)),'o-','linew',1.5,'markersize',6)
plot(hz,spect3(1:length(hz)),'p-','linew',1.5,'markersize',6)
legend({'time average';'coefs average';'power average'})
set(gca,'xlim',[0 25])
xlabel('Frequency (Hz)'), ylabel('Power')

%% Exercise 4

% Rerun the previous simulation using a lot of noise or a little bit of noise (by adjusting the
% multiplication of randn by a large or small number), and larger or small signal (by multiplying
% the signal by a large or small number). Obviously, large signal and low noise is the best scenario,
% but is it better to have large signal but large noise, or small signal and small noise? ”Better” is
% a subjective term, but you can think of it as the visual interpretability of the power spectrum.


% change noise parameter at line 71 and re-run Ex. 2-3 (5 and 500 should be enough)


%% Exercise 5 

% Change the simulation generation from question #2 by creating the three sine waves with random
% phases (between 0 and 2*pi) on each repetition. Then re-run question 3. What do you notice?

srate = 1/1000;
time = 0:srate:2;
pnts = length(time);

freq = [3, 10, 17];
ampl = [5, 7, 4];

for r=1:50

    signal = zeros(1,pnts);

    for i=1:length(ampl)
        signal = signal + ampl(i)*sin(2*pi*freq(i)*time + rand*2*pi);
    end
    
    data(r,:) = signal + 50*randn(size(signal));

end

% plotting
figure(4), clf
subplot(211)
plot(time,signal)
title('Original signal')
subtitle('random phase at each frequency and iteration')

subplot(212)
plot(time,data)
title('Noisy data')


% (1) time average
spect1 = 2*abs(fft( mean(data,1) )/pnts).^2;

% (2) coefficients average
spect2 = 2*abs( mean(fft(data,[],2)/pnts ,1)).^2;

% (3) power average
spect3 = mean( 2*abs( fft(data,[],2)/pnts ).^2,1);

hz = linspace(0,(1/srate)/2,floor(pnts/2)+1);

% plotting
figure(3), clf
plot(hz,spect1(1:length(hz)),'s-','linew',2,'markersize',10)
hold on
plot(hz,spect2(1:length(hz)),'o-','linew',1.5,'markersize',6)
plot(hz,spect3(1:length(hz)),'p-','linew',1.5,'markersize',6)
legend({'time average';'coefs average';'power average'})
set(gca,'xlim',[0 25])
xlabel('Frequency (Hz)'), ylabel('Power')





