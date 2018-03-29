

%signal analysis
[sig, fs] = audioread('D:\1st sem\DSP\DSP Project\ECE6342_morse_code_sound.wav');
len = length(sig);
l = 0 : 1/fs : len/fs - 1/fs;
figure(1);
plot(l, sig);               %time-domain plot in seconds
title('time-domain plot in seconds');
figure(2);
stem(sig);                  %time-domain plot in samples
title('time-domain plot in samples');
sig_fft = abs(fft(sig));
SIG = sig_fft(1:(len/2));
k = 0:fs/(len):fs/2 - fs/(len);
figure(5);
stem(k, SIG);          %freq-domain plot of left side of fft in Hz
title('freq-domain plot of left side of fft in Hz');
j = 0:1/(len/2-1):1;
figure(6);
plot(j, SIG(1:len/2));     %normalized frequency
title('normalized freq of left side of fft in pi radians');
SIG_db = mag2db(SIG);
figure(7);
plot(j, SIG_db(1:len/2));
title('normalized freq of left side of fft in pi radians and magnitude in db');

%designing fir filter conventionally
fn = fs/2;          %Nyquist frequency
fm = 880;          %morsecode frequency
fc = fm/fn;         %normalized morse code frequency
fc1 = fc + 15/fn;           %1st cut off frequency
fc2 = fc - 15/fn;           %2nd cut off frequency
a = fir1(512, [fc2 fc1], 'bandpass');
[H, ~] = freqz(a, floor(len/2));
hold on
plot(j, abs(H), 'r');

n = 3000;
wn = [870 890]/(fs/2);
[~, a] = fir1(n, wn, 'bandpass');

%designing fir filter with Kaiser window
f = [820/fs 825/fs 835/fs 840/fs];
A = [0 1 0];
dev = [0.001 0.05 0.001];
[n,Wn,beta,ftype] = kaiserord(f,A,dev,fs);
b = fir1(n,Wn,'bandpass',kaiser(n+1,beta),'scale');
[H, ~] = freqz(b, floor(len/2));
hold on
plot(j, abs(H), 'r');

%filtering signal
sig_filtered = filtfilt(a, 1, sig);
sig_filtered_fft = abs(fft(sig_filtered));
SIG_filtered = sig_filtered_fft(1:(len/2));
l = 0 : 1/fs : len/fs - 1/fs;
figure(21);
plot(l, sig_filtered);               %time-domain plot in seconds
title('time-domain plot of filtered signal in seconds');
figure(22);
plot(sig_filtered);                  %time-domain plot in samples
title('time-domain plot of filtered signal in samples');
hold on
figure(23);
stem(k, SIG_filtered);                  %time-domain plot in samples
title('fft plot of filtered signal in hz');


%segmenting the signal
seg1 = 13.70;
seg2 = 13.97;
sig_seg = sig(fs*seg1 : fs*seg2);
len_seg = length(sig_seg);
l_seg = seg1 : 1/fs : seg2;
figure(11);
plot(l_seg, sig_seg);               %time-domain plot in seconds
title('time-domain plot in seconds');
sample_seg = fs*seg1 : 1 : fs*seg2;
figure(12);
stem(sample_seg, sig_seg);                  %time-domain plot in samples
title('time-domain plot in samples');
sig_seg_fft = fft(sig_seg);
SIG_seg = abs(sig_seg_fft(1:len_seg/2));
k_seg = 0:fs/(len_seg):fs/2 - fs/(len_seg);
figure(15);
stem(k_seg, SIG_seg);          %freq-domain plot of left side of fft in Hz
title('freq-domain plot of left side of fft in Hz');
j_seg = 0:1/(len_seg/2-1):1;
figure(16);
stem(j_seg, SIG_seg);     %normalized frequency
title('normalized freq of left side of fft in pi radians');
SIG_seg_db = mag2db(SIG_seg);
figure(17);
stem(j_seg, SIG_seg_db(1:len_seg/2));
title('normalized freq of left side of fft in pi radians and magnitude in db');

%plotting filtered segment signal
[H, w] = freqz(a, floor(len/2));
hold on
plot(j, abs(H), 'r');
l_seg = seg1 : 1/fs : seg2;
figure(18);
stem(l_seg, sig_seg_filtered);      %time-domain plot in seconds
k_seg = 0:fs/(len_seg):fs/2 - fs/(len_seg);
figure(19);
stem(k_seg, SIG_seg_filtered);          %freq-domain plot of left side of fft in Hz
title('freq-domain plot of left side of fft in Hz');
j_seg = 0:1/(len_seg/2-1):1;
figure(20);
stem(j_seg, SIG_seg_filtered);     %normalized frequency
title('normalized freq of left side of fft in pi radians');

%butterworth
j = 0:1/(len/2-1):1;
figure(8);
plot(j, SIG_half(1:len/2));     %normalized frequency
SIG_db = mag2db(SIG_half);
hold on
plot(j, SIG_db(1:len/2), 'g');
wp = 1200/fs;
ws = 1600/fs;
Rp = 0.008681;
Rs = 80;
[N, Wn] = buttord(wp,ws,Rp,Rs);
[b, a] = butter(N, Wn, 'low');
H = freqz(b,a, floor(len/2));
hold on
plot(j, abs(H), 'r');
sig_filtered = filter(b, a, sig);
figure(9);
plot(l, sig_filtered);


