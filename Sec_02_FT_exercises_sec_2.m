%% Exercises from section 2

clear, clc

%% Exercise 2

% Write code to generate and plot two-second signals that comprise a sum of the following sine wave
% parameters. Test various sampling rates, ranging from 1 Hz to 1000 Hz, to determine what—if
% any—effect that has on the plots.

% time and srate
srate = 500;
time  = 0:1/srate:2-1/srate;

% define parameters
freqs_1 = [2, 4.2];
phase_1 = [0, pi*0.75];
ampls_1 = [1, 1.7];

freqs_2 = [200, 402, 3.2];
phase_2 = [pi,  0,   1];
ampls_2 = [100, 10, 50];

% define the sinewaves
wave_1 = zeros(size(time));
for i = 1:length(freqs_1)
    wave = ampls_1(i) * sin( 2*pi*freqs_1(i)*time + phase_1(i));
    wave_1 = wave_1 + wave;
end

wave_2 = zeros(size(time));
for i = 1:length(freqs_2)
    wave = ampls_2(i) * sin( 2*pi*freqs_2(i)*time + phase_2(i));
    wave_2 = wave_2 + wave;
end

% plot to check
figure(1), clf
subplot(211)
plot(time,wave_1,'linew',1)
xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
title('Signal 1')
subtitle('sampling rate: 500Hz')

subplot(212)
plot(time,wave_2,'linew',1)
xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
title('Signal 2')
subtitle('sampling rate: 500Hz')

% now try to vary the sampling rate 

srate2 = 1:200:1001;

for j = 1:length(srate2)
    time2  = 0:1/srate2(j):2-1/srate2(j);
    wave_3 = zeros(length(srate2),length(time2));
    for i = 1:length(freqs_1)
        wave = ampls_1(i) * sin( 2*pi*freqs_1(i)*time2 + phase_1(i));
        wave_3 = wave_3 + wave;
    end
    figure(2)
    subplot(2,3,j)
    plot(time2,wave_3,'linew',1)
    xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
    title('Signal 1')
    subtitle(['sampling rate: ' num2str(srate2(j)) 'Hz'])
end

for j = 1:length(srate2)
    time2  = 0:1/srate2(j):2-1/srate2(j);
    wave_4 = zeros(length(srate2),length(time2));
    for i = 1:length(freqs_2)
        wave = ampls_2(i) * sin( 2*pi*freqs_2(i)*time2 + phase_2(i));
        wave_4 = wave_4 + wave;
    end
    figure(3)
    subplot(2,3,j)
    plot(time2,wave_4,'linew',1)
    xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
    title('Signal 1')
    subtitle(['sampling rate: ' num2str(srate2(j)) 'Hz'])
end

%% Exercise 3

% In theory, sine waves have an equal amount of area above and below the Y=0 axis, meaning the
% average of a sine wave over all time points should be 0. Using a sampling rate of 1 kHz, generate
% a 2-second sine wave at 6 Hz and average over all values. It won’t be exactly zero because of
% computer rounding errors, but anything less than 1e-10 you can consider to be zero. 
% Now repeat with a sine wave of 6.1 Hz. Is the answer still zero? Try again with 6.5 Hz. To
% understand why this happens, plot the sine waves.

clear, clc

% time and srate
srate = 1000;
time  = 0:1/srate:2-1/srate;

% parameters
freq = [6, 6.1, 6.5];
amp  = 1;
phs  = 0;

% wave
wave_1 = amp * sin(2*pi*freq(1)*time + phs);
wave_2 = amp * sin(2*pi*freq(2)*time + phs);
wave_3 = amp * sin(2*pi*freq(3)*time + phs);

% take the average of al the value
average1 = mean(wave_1);
average2 = mean(wave_2);
average3 = mean(wave_3);

disp('Signal average:')
disp(average1)
disp('Signal average:')
disp(average2)
disp('Signal average:')
disp(average3)

% plot for clarity
figure(5), clf
plot(time,wave_1,'linew',1.5)
hold on
plot(time,wave_2,'linew',1.5)
plot(time,wave_3,'linew',1.5)
xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
title('Signal')
subtitle('sampling rate: 1000Hz')
legend('6 Hz','6.1 Hz','6.5 Hz')

%% Exercise 4

% Adapt the previous exercise to create and average sine waves in a loop with frequencies varying
% from 1 to 50 Hz in 90 linearly spaced steps. Then make a plot of frequency by average amplitude
% value.

clear, clc

% time and srate
srate = 1000;
time  = 0:1/srate:2-1/srate;

% parameters
freq = 1:(50/91):50;
amp  = 1;
phs  = 0;

% wave
wave    = zeros(size(time));
average = zeros(size(freq));
for i = 1:length(freq)
    wave = amp * sin(2*pi*freq(i)*time + amp);
    average(i) = mean(wave);
end

figure(6)
plot(freq,abs(average),'s--')
set(gca,'ylim',[-.005 .05])
xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
title('Frequency by Average')

%% Exercise 6 

% Draw the following complex numbers on a complex plane. Then extract their magnitude (distance
% to origin) and phase. You can use a calculator to do the inverse tangent (or brush up on your
% trig identities!); the important thing is to write down the equation correctly. Remember that you
% want to square the imaginary component, not the imaginary number itself (thus, (2i)^2
% is treated as 2^2 and does not become -4).

clear, clc

% the formulas (given x a complex number):
%   > |x| = m  =  sqrt( real(x)^2 + imag(x)^2 )
%   > theta(x) = atan2( imag(x)   / real(x)   )

a = complex(2,3);
b = complex(0,1);
c = complex(-3,4);
d = complex(-1,-1);

% compute mag and ang
mag_a = sqrt( real(a)^2 + imag(a)^2 );  % also abs()
ang_a = atan2( imag(a),real(a) );     % also angle()

mag_b = sqrt( real(b)^2 + imag(b)^2 );
ang_b = atan2( imag(b),real(b) );

mag_c = sqrt( real(c)^2 + imag(c)^2 );
ang_c = atan2( imag(c),real(c) );

mag_d = sqrt( real(d)^2 + imag(d)^2 );
ang_d = atan2( imag(d),real(d) );

disp([ 'Mag of a is ' num2str(mag_a) ' and Ang of a is ' num2str(ang_a) '.' ])
disp([ 'Mag of b is ' num2str(mag_b) ' and Ang of b is ' num2str(ang_b) '.' ])
disp([ 'Mag of c is ' num2str(mag_c) ' and Ang of c is ' num2str(ang_c) '.' ])
disp([ 'Mag of d is ' num2str(mag_d) ' and Ang of d is ' num2str(ang_d) '.' ])

% plotting 
figure(7), clf
plot(real(a),imag(a),'s','markersize',12,'linew',2)
hold on
plot(real(b),imag(b),'s','markersize',12,'linew',2)
plot(real(c),imag(c),'s','markersize',12,'linew',2)
plot(real(d),imag(d),'s','markersize',12,'linew',2)

set(gca,'xlim',[-5 5],'ylim',[-5 5])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k','linew',2)
plot([0 0],get(gca,'ylim'),'k','linew',2)
xlabel('Real axis')
ylabel('Imaginary axis')
title('Various complex numbers on the complex plane')
legend(['a: ' num2str(a)],['b: ' num2str(b)],['c: ' num2str(c)],['d: ' num2str(d)])

%% Exercise 7

% Compute the dot product between the following pairs of vectors. Draw the resulting complex dot
% product in complex plane. Finally, compute the magnitude of the dot product.

clear, clc

a = dot( [1,-3], [1,1i] );
b = dot( [4,1], [-2,-2i] );
c = dot( [7i,1], [2,3] );
d = dot( [4,4], [1i,1i] );

disp([ 'Dot product a equals: ' num2str(a)])
disp([ 'Dot product a equals: ' num2str(b)])
disp([ 'Dot product a equals: ' num2str(c)])
disp([ 'Dot product a equals: ' num2str(d)])

% compute mag and ang
mag_a = sqrt( real(a)^2 + imag(a)^2 );  % also abs()
ang_a = atan2( imag(a),real(a) );     % also angle()

mag_b = sqrt( real(b)^2 + imag(b)^2 );
ang_b = atan2( imag(b),real(b) );

mag_c = sqrt( real(c)^2 + imag(c)^2 );
ang_c = atan2( imag(c),real(c) );

mag_d = sqrt( real(d)^2 + imag(d)^2 );
ang_d = atan2( imag(d),real(d) ); 

disp([ 'Mag of a is ' num2str(mag_a) ' and Ang of a is ' num2str(ang_a) '.' ])
disp([ 'Mag of b is ' num2str(mag_b) ' and Ang of b is ' num2str(ang_b) '.' ])
disp([ 'Mag of c is ' num2str(mag_c) ' and Ang of c is ' num2str(ang_c) '.' ])
disp([ 'Mag of d is ' num2str(mag_d) ' and Ang of d is ' num2str(ang_d) '.' ])

% plotting
figure(8), clf
plot(real(a),imag(a),'s','markersize',12,'linew',2)
hold on
plot(real(b),imag(b),'s','markersize',12,'linew',2)
plot(real(c),imag(c),'s','markersize',12,'linew',2)
plot(real(d),imag(d),'s','markersize',12,'linew',2)

set(gca,'xlim',[-15 15],'ylim',[-15 15])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k','linew',2)
plot([0 0],get(gca,'ylim'),'k','linew',2)
xlabel('Real axis')
ylabel('Imaginary axis')
title('Various complex numbers on the complex plane')
legend(['a: ' num2str(a)],['b: ' num2str(b)],['c: ' num2str(c)],['d: ' num2str(d)])



