%% Exercises from section 8

clear, clc

%% Exercise 1

% The two MATLAB examples of side gradients used square-shaped gradients. Adjust the code to
% make the sine gradient a circle instead of a square. Does that impact the resulting amplitude
% spectrum?

% specify vector of sine frequencies
sinefreq = linspace(.0001,.2,50);  % arbitrary units

% leave this fixed for now
sinephas = pi/2;

% sine wave initializations
lims  = [-91 91];
[x,y] = ndgrid(lims(1):lims(2),lims(1):lims(2));
xp    = x*cos(sinephas) + y*sin(sinephas);
clim  = [0 1000];

% setup plot
figure(1), clf
subplot(121)
imageh = imagesc(rand(length(x)));
axis square, axis off, axis xy
title('Space domain')

subplot(222)
amph = imagesc(rand(length(x)));
axis square, axis off, axis xy
set(gca,'xlim',[lims(2)-30 lims(2)+30],'ylim',[lims(2)-30 lims(2)+30],'clim',clim)
title('Amplitude spectrum')

subplot(224)
phaseh = imagesc(rand(length(x)));
axis square, axis off, axis xy
set(gca,'xlim',[lims(2)-30 lims(2)+30],'ylim',[lims(2)-30 lims(2)+30])
title('Phase spectrum')

width = 100;
gaus2d = exp(-(x.^2 + y.^2) ./ (2*width^2));

% now loop over phases
for si=1:length(sinefreq)
    
    % compute sine wave
    img = sin( 2*pi*sinefreq(si)*xp );
    img(gaus2d<.9) = 0;
    
    % 2D FFT and extract power and phase spectra
    imgX  = fftshift(fft2(img));
    
    powr2 = abs(imgX);
    phas2 = angle(imgX);
    
    % update plots
    set(imageh,'CData',img);
    set(amph  ,'CData',powr2);
    set(phaseh,'CData',phas2);
    
    pause(.2)
end

%% Exercise 2

% I showed an example in MATLAB of a solid ball moving around on a plane. The amplitude
% spectrum was a dot in the center. Integrate the sine gradient and ball examples by replacing the
% ball with moving gradients. Then try moving the ball around in different regions of the plane to
% see the effects on the amplitude spectrum.
% It’s also interesting to compare different spatial frequencies. For example, try .01 or .05. For .01,
% the amplitude at low frequencies (center of the amplitude spectrum plane) goes from red to blue
% (decreased amplitude). What’s happening in the image that causes this decrease in amplitude?

lims  = [-91 91];
[x,y] = ndgrid(lims(1):lims(2),lims(1):lims(2));
width = 50;   % width of gaussian

centlocs = linspace(-80,80,50);

% setup plot
figure(3), clf
subplot(121)
imageh = imagesc(rand(length(x)));
axis square, axis off, axis xy
title('Space domain')
set(gca,'clim',[-1 1])

subplot(222)
amph = imagesc(rand(length(x)));
axis square, axis off, axis xy
title('Amplitude spectrum')
set(gca,'clim',[0 200])

subplot(224)
phaseh = imagesc(rand(length(x)));
set(gca,'xlim',[lims(2)-30 lims(2)+30],'ylim',[lims(2)-30 lims(2)+30])
axis square, axis off, axis xy
title('Phase spectrum')

sinefreq = .01;
sinegrad = sin( 2*pi*sinefreq*x );

% now loop over locations (center locations)
for si=1:length(centlocs)
    
    mx = x-centlocs(si);
    my = y-20;
    
    gaus2d = exp(-(mx.^2 + my.^2) ./ (2*width^2));
    img = sinegrad;
    img(gaus2d<.9) = 0;
    
    % 2D FFT and extract power and phase spectra
    imgX  = fftshift(fft2(img));
    powr2 = abs(imgX);
    phas2 = angle(imgX);
    
    % update plots
    set(imageh,'CData',img);
    set(amph  ,'CData',powr2);
    set(phaseh,'CData',phas2);
    
    pause(.2)
end

% for low frequency it happens that the entire 'ball' encompasses less than
% a cycle per 'frame', so that the center of the amplitude spectrum
% oscillates between high anl low values



