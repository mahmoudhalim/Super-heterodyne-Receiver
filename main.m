clc;

addpath('.\ audio');
%% Transmiter
% reading messages files
[message1, Fs] = audioread("Short_BBCArabic2.wav");
[message2, ~] = audioread("Short_FM9090.wav");
% converting message signals to monophonic
message1 = (message1(:,1) + message1(:,2))';
message2 = (message2(:,1) + message2(:,2))';
% extending messages to have the same length
max_length = max(length(message1), length(message2)) + 1;
t = 0:(1/Fs):(max_length-1)/Fs;
message1(1,end+(numel(t)-numel(message1)))=0;
message2(1,end+(numel(t)-numel(message2)))=0;
% Upsampling the audio signals
resampling_factor = 20;
Fs = resampling_factor*Fs;
t = 0:(1/Fs):(max_length*resampling_factor-1)/Fs;
message1_rs = interp(message1, resampling_factor);
message2_rs = interp(message2, resampling_factor);
% modulating the audio signals
Fc = 100000; % carrier frequency
deltaF = 55000; % differnce between carriers
carrier1 = cos(2*pi*Fc*t);
carrier2 =  cos(2*pi*(Fc+deltaF)*t);
modulated_message1 = carrier1 .* message1_rs;
modulated_message2 = carrier2 .* message2_rs;

transmitted_signal = modulated_message1 + modulated_message2;
plotSpectrum(transmitted_signal, Fs, "Magnitude","Transmitter Output");

%% RF Stage
Nc = 0; % desired channel
A_stop1 = 60;		% Attenuation in dB
F_stop1 = (Fc + Nc*deltaF) - 20000;		% Edge of the stopband
F_pass1 = (Fc + Nc*deltaF) - 10000;     % Edge of the passband
F_pass2 = (Fc + Nc*deltaF) + 10000;     % Closing edge of the passband
F_stop2 = (Fc + Nc*deltaF) + 20000;     % Edge of the second stopband
A_stop2 = 60;		% Attenuation in the second stopband in dB
A_pass = 0.1;		% Amount of ripple allowed in the passband in dB

RFbandpassSpecs = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    F_stop1,F_pass1,F_pass2, F_stop2, ...
    A_stop1, A_pass, A_stop2, Fs);

RFbandpassFilter = design(RFbandpassSpecs,"equiripple");
% fvtool(bandpassFilter)
RF_output = filter(RFbandpassFilter, transmitted_signal);
plotSpectrum(RF_output, Fs,"Magnitude","RF Stage Output");

%% Mixer
IF = 27500 ; % Intermidiate frequency
Flo = Fc+Nc*deltaF+ IF; % LO frequency
IF_carrier =cos(2*pi*Flo*t);
Mixer_output = RF_output .* IF_carrier;
plotSpectrum(Mixer_output, Fs, "Magnitude", "RF Mixer Output");

%% IF Stage
A_stop1 = 60;		% Attenuation in dB
F_stop1 = IF - 20000;		% Edge of the stopband
F_pass1 = IF - 10000;     % Edge of the passband
F_pass2 = IF + 10000;     % Closing edge of the passband
F_stop2 = IF + 20000;     % Edge of the second stopband
A_stop2 = 60;		% Attenuation in the second stopband in dB
A_pass = 0.1;		% Amount of ripple allowed in the passband in dB

IFbandpassSpecs = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    F_stop1,F_pass1,F_pass2, F_stop2, ...
    A_stop1, A_pass, A_stop2, Fs);

IFbandpassFilter = design(IFbandpassSpecs,"equiripple");
% fvtool(bandpassFilter)
IF_output = filter(IFbandpassFilter, Mixer_output);
plotSpectrum(IF_output, Fs, "Magnitude", "IF Stage Output");

%% Baseband Detection
BB_carrier = cos(2*pi*IF*t);
BB_message = BB_carrier .* IF_output;
plotSpectrum(BB_message, Fs, "Magnitude", "Output of Mixer");
F_pass = 15000; % Edge of the lowband
F_stop = 20000; % Edge of the stopband
A_stop = 60; % Attenuation in the band

% low pass filter
BBlowpassSpecs = fdesign.lowpass('Fp,Fst,Ap,Ast', ...
    F_pass, F_stop, A_pass, A_stop, Fs);
BBlowpassfilter = design(BBlowpassSpecs, "butter");
% fvtool(BBlowpassfilter);
output_message = filter(BBlowpassfilter, BB_message);
plotSpectrum(output_message, Fs,"Magnitude", "Output of the LPF");

output_message = decimate(output_message, resampling_factor); % down sample the messsage to its original fs
audiowrite("out1.wav",output_message, Fs/10);
%% Helpers;
function [] = plotSpectrum(y, fs, label, Ftitle)
N = length(y);
Y = fft(y,N);
f = (-N/2:N/2-1)*fs/N;
figure;
plot(f,abs(fftshift(Y))/N);
xlabel("Frequency (Hz)");
ylabel (label);
title(Ftitle)

end
