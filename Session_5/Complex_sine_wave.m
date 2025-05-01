srate = 250;
time = 0:1/srate:1;
Freq = 5;

real_sin = sin(2*pi*Freq*time);
complex_sin = exp(1i*2*pi*(Freq-1).*time);

figure(1); clf

subplot(311);
plot(time, real_sin)
title("Real sine wave")

subplot(312);
plot(time, imag(complex_sin))
title("Imaginary complex sine wave") % Sine component

subplot(313);
plot(time, real(complex_sin))
title("Real complex sine wave") %phase shifted (cosine component)

figure(2); clf
plot(time, real(complex_sin))
hold on
plot(time, imag(complex_sin), 'r--')

figure(3); clf
plot3(time, imag(complex_sin), real(complex_sin));

%% How real and imaginary part could unravel amplitude and phase

% you can not extract phase from one point. but with complex sinewaves you
% can easily find the phase by finding tha angle between real and imaginary
% part and abs of both

A1 = 3 * cos(2*pi*Freq*time + 0);
A2 = 2 * cos(2*pi*Freq*time + pi/6);
A3 = 1 * cos(2*pi*Freq*time + pi/3);

fCoef1 = fft(A1)/length(time);
fCoef2 = fft(A2)/length(time);
fCoef3 = fft(A3)/length(time);

hz = linspace(0,srate/2,floor(length(time)/2) + 1);
hz6 = dsearchn(hz',Freq);

figure;
h(1) = polar([0 angle(fCoef1(hz6))], [0 2*abs(fCoef1(hz6))], 'r');
hold on
h(2) = polar([0 angle(fCoef2(hz6))], [0 2*abs(fCoef2(hz6))], 'b');
h(3) = polar([0 angle(fCoef3(hz6))], [0 2*abs(fCoef3(hz6))], 'y');
set(h, 'linewidth', 5);
legend({'signal 1'; 'signal2';'signal3'})



