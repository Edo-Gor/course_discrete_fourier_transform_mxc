%% Exercises from section 4

clear, clc

%% Exercise 1

% Open a new script in MATLAB or Python and program the loop-based discrete inverse Fourier
% transform from scratch. Test it by generating a random vector, taking the forward Fourier
% transform, then the inverse Fourier transform, and then compare the original to the 
% reconstructed signal

% compute the direct FT
srate = 1/1000;
time = 0:srate:3;
pnts = length(time);

freq1 = 5; freq2 = 8;
ampl1 = 2; ampl2 = 3;

wave_signal = ampl1*sin(2*pi*freq1*time) + ampl2*sin(2*pi*freq2*time);
rand_signal = randn(1,pnts);

FourTime = (0:pnts-1) / pnts;

Wcoeffs = zeros(size(time));
Rcoeffs = zeros(size(time));

for i = 1:pnts
    
    csw = exp(-1i*2*pi*(i-1)*FourTime);
    
    Wcoeffs(i) = sum(csw.*wave_signal);
    Rcoeffs(i) = sum(csw.*rand_signal);
    
end

Wamplitude = abs(Wcoeffs/pnts);
Wamplitude(2:end) = 2*Wamplitude(2:end);

Ramplitude = abs(Rcoeffs/pnts);
Ramplitude(2:end) = 2*Ramplitude(2:end);

hz = linspace(0, (1/srate)/2, floor(pnts/2)+1);


% plotting wave signal
figure(1), clf
subplot(211)
plot(time,wave_signal,'k','linew',2)
xlabel('Time (ms)')
ylabel('Amplitude')
title('Time domain (composite sinewave)')

subplot(212)
stem(hz,Wamplitude(1:length(hz)),'ks-','linew',2,'markerfacecolor','w','markersize',10)
hold on

fCoefsF = fft(wave_signal)/length(time);
amplsF  = 2*abs(fCoefsF);
stem(hz,amplsF(1:length(hz)),'ro','markerfacecolor','r')

set(gca,'xlim',[0 20],'ylim',[0,3.5])
legend({'Amplitude (manual loop)';'Amplitude (fft() function)'})
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Frequency domain')


% plotting random signal
figure(2), clf
subplot(211)
plot(time,rand_signal,'k','linew',2)
xlabel('Time (ms)')
ylabel('Amplitude')
title('Time domain (random signal)')

subplot(212)
stem(hz,Ramplitude(1:length(hz)),'ks-','linew',2,'markerfacecolor','w','markersize',10)
hold on

fCoefsF = fft(rand_signal)/pnts;
amplsF  = 2*abs(fCoefsF);
stem(hz,amplsF(1:length(hz)),'ro','markerfacecolor','r')

set(gca,'xlim',[0 20],'ylim',[0,0.2])
legend({'Amplitude (manual loop)';'Amplitude (fft() function)'})
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Frequency domain')

% compute the inverse FT
rec_wave_signal = zeros(size(wave_signal));
rec_rand_signal = zeros(size(rand_signal));

for j = 1:pnts
    
    csw = exp(1i*2*pi*(j-1)*FourTime);
    
    w_csw = csw * Wcoeffs(j);
    r_csw = csw * Rcoeffs(j);
    
    rec_wave_signal = rec_wave_signal + w_csw;
    rec_rand_signal = rec_rand_signal + r_csw;
    
end

rec_wave_signal = rec_wave_signal / pnts;
rec_rand_signal = rec_rand_signal / pnts;

% plotting the reconstructed signals
figure(3), clf
subplot(211)
plot(time,wave_signal)
hold on
plot(time,real(rec_wave_signal),'ro')
legend({'original';'reconstructed'})
xlabel('Time (s)')

subplot(212)
plot(time,rand_signal)
hold on
plot(time,real(rec_rand_signal),'ro')
legend({'original';'reconstructed'})
xlabel('Time (s)')

%% Exercise 2

clear, clc

% In a previous section, I mentioned two normalizations to apply to the output of the Fourier
% transform when interpreting amplitude or power. Do you still need to normalize when applying
% the inverse Fourier transform? To find out, generate a signal, take its Fourier transform, and then
% apply the inverse Fourier transform. Compute the inverse twice: Once without normalizing the
% Fourier coefficients before computing the inverse, and once with normalizations before computing
% the inverse.

srate = 1/1000;
time = 0:srate:3;
pnts = length(time);
freq1 = 5; freq2 = 8;
ampl1 = 2; ampl2 = 3;

signal = ampl1*sin(2*pi*freq1*time) + ampl2*sin(2*pi*freq2*time);

fCoefsF1 = fft(signal)/pnts;
amplsF1  = 2*abs(fCoefsF1);

fCoefsF2 = fft(signal);

rec_norm_signal = ifft(amplsF1);
rec_nnor_signal = ifft(fCoefsF2);

% plotting the reconstructed signals
figure(4), clf
subplot(211)
plot(time,signal)
hold on
plot(time,real(rec_norm_signal),'ro')
legend({'original';'reconstructed with normlisation'})
xlabel('Time (s)')

subplot(212)
plot(time,signal)
hold on
plot(time,real(rec_nnor_signal),'ro')
legend({'original';'reconstructed without normalisation'})
xlabel('Time (s)')

%% Exercises 3-5

clear, clc

% Create a 1-second sine wave at 5 Hz using 100 Hz sampling rate. This signal will be used for
% the remaining exercises. Compute the Fourier transform and obtain two copies of the Fourier
% coefficients—one for plotting, and the other for modulating. Then plot the time-domain signal
% in the top panel and the amplitude spectrum on the bottom panel.

srate = 1/100;
time = 0:srate:1-srate;
pnts = length(time);
freq1 = 5; 
ampl1 = 1;

signal = ampl1*sin(2*pi*freq1*time);

fCoefsF1 = fft(signal);
fCoefsF2 = fft(signal);

hz = linspace(0,(1/srate)/2,floor(pnts/2)+1);

% plotting
figure(5), clf
subplot(211)
plot(time,signal,'linew',2)
title('Time domain signal')
xlabel('Time (s)')

subplot(212)
stem(1:pnts,abs(fCoefsF1),'s-','linew',2,'markersize',10,'markerfacecolor','w')
title('Frequency spectrum')
xlabel('Time (s)')

% Change the phase of the 5 Hz component while preserving the amplitude. To do this, you need
% to transform the 5 Hz Fourier coefficient to its Euler notation (ae^{ith}), then change th and then
% replace the 5 Hz coefficient in the extra copy of the Fourier coefficients. (Hint: don’t forget
% about the negative frequencies!) Plot the two time-domain signals (original and reconstructed)
% on top of each other.

new_phase = pi/8;
[~, fidx] = min(abs(hz-freq1));

fCoefsF2(fidx) = abs(fCoefsF1(fidx)) * exp(1i*new_phase);
fCoefsF2(end-fidx+2) = abs(fCoefsF1(end-fidx+2)) * exp(-1i*new_phase);

rec_signal = real(ifft(fCoefsF2));

% plotting
figure(6), clf
subplot(311)
plot(time,signal)
hold on
plot(time,rec_signal,'ro')
title('Positive and negative phase modulation (pi/8)')
legend({'original';'rec. after phase modulation'})
xlabel('Time (s)')
set(gca,'xlim',[0 1],'ylim',[-1.05,1.05])


% Comment out the line where you modulate the negative frequency component, and set the phase
% of +5 Hz to be pi/2. What happens to the reconstructed signal and why? Try it again with the
% phase set to pi/4.

new_phase = pi/2;

fCoefsF2 = fft(signal);
fCoefsF2(fidx) = abs(fCoefsF1(fidx)) * exp(1i*new_phase);
%fCoefsF2(end-fidx+1) = abs(fCoefsF1(end-fidx+1)) * exp(-1i*new_phase);

rec_signal = real(ifft(fCoefsF2));

% plotting
subplot(312)
plot(time,signal)
hold on
plot(time,rec_signal,'ro')
title('Positive phase modulation only (pi/2)')
legend({'original';'rec. after phase modulation'})
xlabel('Time (s)')
set(gca,'xlim',[0 1],'ylim',[-1.05,1.05])


new_phase = pi/4;

fCoefsF2 = fft(signal);
fCoefsF2(fidx) = abs(fCoefsF1(fidx)) * exp(1i*new_phase);
%fCoefsF2(end-fidx+1) = abs(fCoefsF1(end-fidx+1)) * exp(-1i*new_phase);

rec_signal = real(ifft(fCoefsF2));

% plotting
subplot(313)
plot(time,signal)
hold on
plot(time,rec_signal,'ro')
title('Positive phase modulation only (pi/4)')
legend({'original';'rec. after phase modulation'})
xlabel('Time (s)')
set(gca,'xlim',[0 1],'ylim',[-1.05,1.05])

%% Exercises 6-7

% Repeat questions 4 and 5 but simulate the signal for 0.9 seconds instead of 1.0 seconds. How
% are the results different and why?

srate = 1/100;
time = 0:srate:0.9-srate;
pnts = length(time);
freq1 = 5; 
ampl1 = 1;

signal = ampl1*sin(2*pi*freq1*time);

fCoefsF1 = fft(signal);
fCoefsF2 = fft(signal);

hz = linspace(0,(1/srate)/2,floor(pnts/2)+1);

% plotting
figure(7), clf
subplot(211)
plot(time,signal,'linew',2)
title('Time domain signal (incomplete cicle)')
xlabel('Time (s)')

subplot(212)
stem(1:pnts,abs(fCoefsF1),'s-','linew',2,'markersize',10,'markerfacecolor','w')
title('Frequency spectrum')
xlabel('Time (s)')

% change phases

new_phase = pi/8;
[~, fidx] = min(abs(hz-freq1));

fCoefsF2(fidx) = abs(fCoefsF1(fidx)) * exp(1i*new_phase);
fCoefsF2(end-fidx+2) = abs(fCoefsF1(end-fidx+2)) * exp(-1i*new_phase);

rec_signal = real(ifft(fCoefsF2));

% plotting
figure(8), clf
subplot(311)
plot(time,signal)
hold on
plot(time,rec_signal,'ro')
title('Positive and negative phase modulation (pi/8)')
legend({'original (incomplete cicle)';'rec. after phase modulation'})
xlabel('Time (s)')
set(gca,'xlim',[0 1],'ylim',[-1.05,1.05])

new_phase = pi/2;

fCoefsF2 = fft(signal);
fCoefsF2(fidx) = abs(fCoefsF1(fidx)) * exp(1i*new_phase);
%fCoefsF2(end-fidx+1) = abs(fCoefsF1(end-fidx+1)) * exp(-1i*new_phase);

rec_signal = real(ifft(fCoefsF2));

% plotting
subplot(312)
plot(time,signal)
hold on
plot(time,rec_signal,'ro')
title('Positive phase modulation only (pi/2)')
legend({'original (incomplete cicle)';'rec. after phase modulation'})
xlabel('Time (s)')
set(gca,'xlim',[0 1],'ylim',[-1.05,1.05])

new_phase = pi/4;

fCoefsF2 = fft(signal);
fCoefsF2(fidx) = abs(fCoefsF1(fidx)) * exp(1i*new_phase);
%fCoefsF2(end-fidx+1) = abs(fCoefsF1(end-fidx+1)) * exp(-1i*new_phase);

rec_signal = real(ifft(fCoefsF2));

% plotting
subplot(313)
plot(time,signal)
hold on
plot(time,rec_signal,'ro')
title('Positive phase modulation only (pi/4)')
legend({'original (incomplete cicle)';'rec. after phase modulation'})
xlabel('Time (s)')
set(gca,'xlim',[0 1],'ylim',[-1.05,1.05])

% To understand these effects better, make a plot that shows the Fourier coefficients in the complex
% plane. Show the positive and negative coefficients for the original and modified spectra.

fCoefsF1 = fft(signal);
fCoefsF2 = fft(signal);

new_phase = pi/8;
[~, fidx] = min(abs(hz-freq1));
fCoefsF2(fidx) = abs(fCoefsF1(fidx)) * exp(1i*new_phase);
fCoefsF2(end-fidx+2) = abs(fCoefsF1(end-fidx+2)) * exp(-1i*new_phase);

fCoefsF1 = fCoefsF1/length(time);
fCoefsF2 = fCoefsF2/length(time);

figure(9), clf
plot(real(fCoefsF1(fidx)),imag(fCoefsF1(fidx)),'ro','markerfacecolor','r','markersize',12)
hold on
plot(real(fCoefsF1(end-fidx+2)),imag(fCoefsF1(end-fidx+2)),'rs','markerfacecolor','r','markersize',12)
plot(real(fCoefsF2(fidx)),imag(fCoefsF2(fidx)),'ko','markerfacecolor','k','markersize',12)
hold on
plot(real(fCoefsF2(end-fidx+2)),imag(fCoefsF2(end-fidx+2)),'ks','markerfacecolor','k','markersize',12)

axis([-1 1 -1 1]), grid on
plot(get(gca,'xlim'),[0 0],'k')
plot([0 0],get(gca,'xlim'),'k')
legend({'Orig +ve';'Orig -ve';'Mod +ve';'Mod -ve'})
axis square
xlabel('Real'), ylabel('Imag')
title('Fourier coefficients in the complex plane')
