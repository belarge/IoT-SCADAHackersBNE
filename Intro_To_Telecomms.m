% This is the exercise files to be used for the IoT and SCADA Hackers session presented April 2020

pkg load signal

t = 0:0.001:1; %define 1 seconds at 0.001 second increases this is effectively 1000hz sampling rate
t2 = 0:0.1:1;
f = 40 %define the frequency at 40 Hz
freq = -500:1:500; %crates a vector of frequency values - Nyquist sampling of 1000hz/2 -> Vector of -500 to 500 Hz 


%Section 1 - The Time Domain (slide 14)
%*****************************************************

%plot a simple cosine function
y = cos (2*pi*f*t);
plot(t(1:100),y(1:100)) #plotting only the first 100 values of signal y

%Defining the chart info
title("Plot of a Sinusoid")
xlabel("time")
axis([0 0.1])

y2 = cos(2*pi*f*t) + 0.6*cos(2*pi*f*t*3) + 0.2*sin(2*pi*f*t*5);
plot(t(1:100),y2(1:100)) #plotting only the first 100 values of signal y


%Defining the chart info
title("Plot of a y_2")
xlabel("time")
axis([0 0.1])

%hold example
clf
hold
plot(cos(2*pi*f*t)(1:100))
plot(0.6*cos(2*pi*f*t*3)(1:100))
plot(0.2*sin(2*pi*f*t*5)(1:100))
plot(y2(1:100),"o-")

%build up more complicated waveforms
%plot a square wave
square_wave = square([1:100])
plot(square_wave(1:30))
title("Plot of a square wave")
xlabel("time")
axis([0 30 0 1.5])
%plot a sawtooth wave
clf
sawtooth_wave=[sawtooth(t2),sawtooth(t2),sawtooth(t2),sawtooth(t2),sawtooth(t2),sawtooth(t2),sawtooth(t2)]
plot(sawtooth_wave)
title("Plot of a sawtooth wave")
xlabel("time")



%Section 2 - The frequency domain (Slide 21)
%*****************************************************


%complex numbers revision in the tool 

complex_example = 2+3j

disp ("Complex Number revision for ... ");
complex_example
disp ("Real Component");
real(complex_example)

disp ("Immaginary Component");
imag(complex_example)

disp ("|complex_example| - The Magnitude of the imaginary number")
abs(complex_example)

disp ("The Phase of the imaginary number")
angle(complex_example)


%FFT

plot(y(1:100))

clf
Y_FFT=fftshift(fft(y));
Y2_FFT = fftshift(fft(y2));

plot(freq,abs(fftshift(fft(y))))
title("FFT of y")
xlable("Frequency - Hz") 

#weird example
plot(y2(1:100))

plot(fft(y2))


plot(freq,abs(fftshift(fft(y2))))
title("Plot of FFT of y_2")
xlabel("Frequency - Hz")

% no FFT shift
plot(abs(fft(y2)))

plot(freq,fftshift(abs(fft(y2))))
title("Plot of FFT of y_2")
xlabel("Frequency - Hz")

clf
subplot(1,2,1)
plot(freq,fftshift(abs(fft(y))))
title("Plot of FFT of y")
xlabel("Frequency - Hz")

subplot(1,2,2)
plot(freq,fftshift(abs(fft(y2))))
title("Plot of FFT of y_2")
xlabel("Frequency - Hz")


#Section 3 - Filtering

%types of filters - LPF, HPF, BPF
Low_Pass_Filter = [zeros(1,330),ones(1,340),zeros(1,331)]; %let's build a low pass filter
Low_Pass_Filter1 = [zeros(1,309),0:0.05:1,ones(1,340),1:-0.05:0,zeros(1,310)]; %let's build a low pass filter
High_Pass_Filter = 1-Low_Pass_Filter;
Band_Pass_Filter = [zeros(1,150),ones(1,250),zeros(1,200),ones(1,250),zeros(1,151)];

#examples of filters 
clf
subplot(3,1,1)
plot(freq,Low_Pass_Filter)
title("Plot of a Low Pass Filter")
xlabel("Frequency - Hz")
axis([-500 500 0 1.5])

subplot(3,1,2)
plot(freq,High_Pass_Filter)
title("Plot of a High Pass Filter")
xlabel("Frequency - Hz")
axis([-500 500 0 1.5])

subplot(3,1,3)
plot(freq,Band_Pass_Filter)
title("Plot of a Band Pass Filter")
xlabel("Frequency - Hz")
axis([-500 500 0 1.5])

#example of the BPF
clf
plot(freq,abs(fftshift(fft(y2))))
hold
plot(freq,Band_Pass_Filter.*550) %dirty hack
axis([-500 500 0 600])

%to do 


#Section 4 - Modulation (slide 31)

%putting the concepts together let's talk about modulation, and first up Amplitude Modulation 
carrier = cos(2*pi*300*t);

%let's modulate y2 onto carrier
modulated_signal = 0.5*y2.*carrier;
clf
subplot(2,2,1)
plot(freq,abs(fftshift(fft(carrier))))
title("Frequency Domain of carrier")
xlabel("Frequency (Hz)") 
subplot(2,2,2)
plot(freq,abs(fftshift(fft(y2))))
title("Frequency Domain of Signal - y_2")
xlabel("Frequency (Hz)") 
subplot(2,1,2)
plot(freq,abs(fftshift(fft(modulated_signal))))
title("Frequency Domain of Modulated Signal")
xlabel("Frequency (Hz)") 


