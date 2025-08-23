%% Exercises from section 3

clear, clc

%% Exercise 2

% By now, you’ve read the description of how the Fourier transform works and have seen it in code.
% Open a new script in MATLAB and program the loop-based discrete Fourier transform from scratch. 
% Test it using a simple signal, like a sine wave, where you can confirm the accuracy of your code.

srate = 1/1000;
time = 0:srate:2;
pnts = length(time);
freq = 12; freq2 = 5;
ampl = 2;  ampl2 = 1;

signal = ampl*sin(2*pi*freq*time) + ampl2*sin(2*pi*freq2*time);

coeffs = zeros(size(signal));
FourierTime = (0:pnts-1)/pnts;

for i = 1:pnts
    
    csw = exp(-1i*2*pi*(i-1)*FourierTime);
    coeffs(i) = sum(signal.*csw) / pnts;
    
end

signal_amplitude = abs(coeffs);
signal_amplitude(2:end) = 2*signal_amplitude(2:end);

hz = linspace(0, (1/srate)/2, floor(pnts/2)+1);

% plotting
figure(1), clf
subplot(211)
plot(time,signal,'k','linew',2)
xlabel('Time (ms)')
ylabel('Amplitude')
title('Time domain')

subplot(212)
stem(hz,signal_amplitude(1:length(hz)),'ks-','linew',2,'markerfacecolor','w','markersize',10)
hold on

fCoefsF = fft(signal)/length(time);
amplsF  = 2*abs(fCoefsF);
stem(hz,amplsF(1:length(hz)),'ro','markerfacecolor','r')

set(gca,'xlim',[0 20],'ylim',[0,2.5])
legend({'Amplitude (manual loop)';'Amplitude (fft() function)'})
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Frequency domain')

%% Exercise 6

% You can compute power either as abs(x)^2 or x*conj(x). I use the former 
% when teaching because of the clear link to geometry. But is one faster than the other? 
% To find out, generate a 10,000-point random signal, compute the Fourier transform, 
% and use these two methods to compute power. Because this is a scientific experiment, 
% you should run the procedure many times on different vectors, and then average 
% the resulting computation times together. Show the resulting average clock times in a bar plot.

time = zeros(2,100);

for i=1:100
    sx = fft(randn(10000,1));
    tic; powr=abs(sx).^2; time(1,i)=toc;
    tic; powr=sx.*conj(sx); time(2,i)=toc;
end

figure(2), clf
bar(mean(time,2)), hold on
errorbar(mean(time,2),std(time,[],2),'o','linew',1.5)
ylabel('Computation time (s)')
set(gca,'xlim',[0 3],'xticklabel',{'abs(sx)^2';'sx.*conj(sx)'})




