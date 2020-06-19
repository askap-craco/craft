% -----------------------------------------------------------------------------
%  (c) Copyright - Commonwealth Scientific and Industrial Research Organisation
%  (CSIRO) - 2016
% 
%  All Rights Reserved.
% 
%  Restricted Use.
% 
%  Copyright protects this code. Except as permitted by the Copyright Act, you
%  may only use the code as expressly permitted under the terms on which the
%  code was licensed to you.
% 
% -----------------------------------------------------------------------------
%    File Name:             ADE_analysis_synthesis_WOLA.m
%    Type:                  Matlab function
%    Contributing authors:  J. Tuthill
%    Created:               Sun 3 April 08:14:48 2016
% 
%    Title:                 Oversampled polyphase analysis/synthesis filterbanks
%    Description:           Matlab function implements an oversampled
%                           polyphase analysis filterbank followed by an
%                           oversampled synthesis filterbank based on the
%                           Weighted-Overlap-Add method.
%
% -----------------------------------------------------------------------------

clear all;

%% Filterbank prototype filters
% analysis filter

%Load the prototype filter for the ADE 6-lane mixed-radix analysis filterbank.
load ADE_R6_OSFIR;
K = length(c); % length of prototype filter
h1 = c;
paths = 1536; % number of polyphase paths (transform length)
taps = K/paths; % length of each polyphase FIR filter
num = 32;
den = 27;
OS = num/den; % oversampling ratio
stride = paths/OS; % stride length through the input data per frame
ovlp = paths-stride; % overlap per frame
fft_shifts = mod((0:num-1)*(paths-stride),paths);
fft_shifts_rev = fliplr(fft_shifts);

h1a=h1.*exp(1i*2*pi*(0:K-1)*1/paths);
h1b=h1.*exp(1i*2*pi*(0:K-1)*2/paths);
h1os=h1.*exp(1i*2*pi*(0:K-1)*OS/paths);

% synthesis filter
h2=h1;

figure(1)
clf
subplot(2,1,1)
plot(h1,'linewidth',2)
hold on
plot(h2,'r','linewidth',2)
hold off
grid on
title(sprintf('Impulse Response, %d-Tap Prototype Nyquist Filter for\n%d-Path Analysis Channelizer, %d-Taps per Path',K, paths, taps))
xlabel('Time Index')
ylabel('Amplitude')

NFFT = 100000;
freq_ax = (-NFFT/2:NFFT/2-1)/NFFT*paths;
subband_0 = fftshift(20*log10(abs(fft(h1/sum(h1),NFFT))));
subband_1 = fftshift(20*log10(abs(fft(h1a/sum(h1),NFFT))));
subband_2 = fftshift(20*log10(abs(fft(h1b/sum(h1),NFFT))));
plot_range = (49800:50500);
subplot(2,1,2)
plot(freq_ax(plot_range),subband_0(plot_range))
hold on
plot(freq_ax(plot_range),subband_1(plot_range),'r')
plot(freq_ax(plot_range),subband_2(plot_range),'g')
hold off
grid on
axis([-2 5 -100 10])
title(sprintf('Frequency Response, Prototype Nyquist Filter response\nfor a %d-Path Analysis Channelizer',paths))
xlabel('Frequency (MHz)')
ylabel('Log Magnitude (dB)')

plot_range = (49700:50300);
subband_0 = fftshift(20*log10(abs(fft(h1/sum(h1),NFFT))));
% image_1p = fftshift(20*log10(abs(fft(h1b/sum(h1),NFFT))));
% image_1m = fftshift(20*log10(abs(fft(conj(h1b)/sum(h1),NFFT))));
image_osp = fftshift(20*log10(abs(fft(h1os/sum(h1),NFFT))));
image_osm = fftshift(20*log10(abs(fft(conj(h1os)/sum(h1),NFFT))));
synth_filt = fftshift(20*log10(abs(fft(h2/sum(h2),NFFT))));
figure(2)
clf
plot(freq_ax(plot_range),subband_0(plot_range),'linewidth',2)
hold on
% plot(freq_ax(plot_range),image_1p(plot_range),'linewidth',2)
% plot(freq_ax(plot_range),image_1m(plot_range),'linewidth',2)
plot(freq_ax(plot_range),image_osp(plot_range),'g','linewidth',2)
plot(freq_ax(plot_range),image_osm(plot_range),'g','linewidth',2)
plot(freq_ax(plot_range),synth_filt(plot_range),'r','linewidth',2)
hold off
grid on
% axis([-3 3 -100 10])
title(sprintf('Frequency Response, Prototype Nyquist Filter for:\n%d-Path Analysis Channelizer (Blue) and\n image bands (green)\nfor %d-Path Synthesis Channelizer (Red)',paths,paths))
xlabel('Frequency (MHz)')
ylabel('Log Magnitude (dB)')

%% Filterbank input signal
fs = 1536e6; %sampling frequency
subband_id = 201;
K = length(c);
N = 100*K; % length of input sequence
N = floor(N/stride)*stride;

t = (0:N-1)/fs;
a_sig_1 = 1; % tone 1 amplitude
% f_sig_1 = 500*1e6; % test 1 tone frequency
f_sig_1 = (subband_id-1)*1e6; % test 1 tone frequency
a_sig_2 = 0.5; % tone 2 amplitude
% f_sig_2 = 500.2e6; % test 2 tone frequency
f_sig_2 = (subband_id-1)*1e6 + 200e3; % test 2 tone frequency

snr = 10; % signal-to-noise ratio

a_noise = 10^((20*log10(a_sig_1/sqrt(2)) - snr)/10);
in_noise = sqrt(a_noise)*randn(1,N);

pos = 0;
in = [zeros(1,pos) 1 zeros(1,N-pos-1)];
in = a_sig_1*cos(2*pi*f_sig_1*t) + in_noise;
% in = a_sig_1*cos(2*pi*f_sig_1*t);
in = in + a_sig_2*cos(2*pi*f_sig_2*t);
% in = a_sig_1*cos(2*pi*f_sig_2*t) + in_noise;

load BETA_bf_filt; % load a basic bandpass filter (models the analogue anti-aliasing filter)
in_filt = filter(Num,1,in);
x = in_filt;

in_spect = 20*log10(abs(fft(in_filt)));
in_spect = in_spect - max(in_spect);
freq_ax = (0:length(in_filt)/2-1)*fs/length(in_filt);
figure(3);
clf
plot(freq_ax/1e6,in_spect(1:length(freq_ax)))
xlabel('Frequency (MHz)')
ylabel('Magnitude (dB)')

%% Analysis-synthesis filterbanks

N = length(x);
num_frames = floor((N-K)/stride);

hh1 = reshape(fliplr(h1),paths,taps); % reshape the analysis filter into the polyphase structure
v1 = zeros(1,paths)';
af_in_vec = zeros(1,paths*taps);
af_in_vec(1:K-stride) = fliplr(x(1:K-stride));
offset = K-stride;
reg1 = zeros(paths,taps);
af_out = zeros(paths,num_frames);

% Analysis filterbank section
for frame = 0:num_frames-1
    
    af_in_vec(stride+1:K) = af_in_vec(1:K-stride);
    af_in_vec(1:stride) = fliplr(x(offset+frame*stride+1:offset+frame*stride+stride));
    reg1 = reshape(af_in_vec,paths,taps);
    
    for k=1:paths
        v1(k)=reg1(k,:)*hh1(k,:)';
    end
    v1 = flipud(v1);
    
    shift_indx = fft_shifts(mod(frame,length(fft_shifts))+1);
    v1 = [v1(shift_indx+1:end);v1(1:shift_indx)];
    
    v2 = paths*fft(v1);
    af_out(:,frame+1) = v2;
    
end

% Synthesis filterbank section
sf_in = af_out.';
s = size(sf_in);
num_frames = s(1);%number of FFT blocks
reg = zeros(1,paths*taps); % synthesis filterbank register
sf_state = zeros(1,K+stride*(num_frames-1)); % output data drray
sf_out = [];
for frame = 0:num_frames-1
    
    % Weighted Overlap-Add Method
    % Periodically-extend the input frame to fill the filter ragister
    reg = zeros(1,paths*taps);
    in_frame = sf_in(frame+1,:); %input frame is 'paths' long
%     in_frame(subband_id-20:subband_id+20) = 0;
    temp = ifft(in_frame);
    % Oversampling correction circular shift
    shift_indx = fft_shifts_rev(mod(frame,length(fft_shifts))+1);
    temp = [temp(shift_indx+1:end) temp(1:shift_indx)];
    for lp1 = 0:taps-1
        reg(1+lp1*paths:(lp1+1)*paths) = temp;
    end
    reg = reg.*h2;
    range = (1+frame*stride):(K + frame*stride);
    sf_state(range) = sf_state(range)+reg;
    if frame >  0
        sf_out = [sf_out sf_state(1+(frame-1)*stride:frame*stride)]; %output frame is 'stride' long
    end
end

out_nr = os_pfb(x, h1, paths, stride, fft_shifts);

subband_time_series = af_out(subband_id,:);
NFFT = length(subband_time_series);
mag_1 = 20*log10(abs(fftshift(fft(subband_time_series))));
mag_1 = mag_1 - max(mag_1);
mag_2 = 20*log10(abs(fftshift(fft(out_nr(:,subband_id))))); % extend the sub-band time-series output from the os_pfb() function to match the time-series length from the WOLA synthesis filter
mag_2 = mag_2 - max(mag_2);
freq_ax = (-NFFT/2:NFFT/2-1)/NFFT*1e6*OS;
figure(4)
clf
plot(freq_ax/1e6,mag_1)
hold on
plot(freq_ax/1e6,mag_2,'r')
hold off
xlabel('Frequency (MHz)')
ylabel('Magnitude (dB)')
title(sprintf('Sub-band power spectrum\nSub-band number = %d',subband_id))

mag_1 = 20*log10(abs(fft(sf_out)));

mag_1 = mag_1 - max(mag_1);
NFFT = length(sf_out);
fax = (0:NFFT/2-1)*fs/NFFT;
figure(5)
clf
plot(fax/1e6,mag_1(1:length(mag_1)/2))
xlabel('Frequency (MHz)')
ylabel('Magnitude (dB)')
