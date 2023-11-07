clear;

[s1, fs1] = audioread("BASIC5000_0269.wav");
[s2, fs2] = audioread("TRAVEL1000_0282.wav");

assert(fs1 == fs2);
fs = fs1;



max_length = max(length(s1), length(s2));


s1 = [s1 ;zeros(max_length-length(s1), 1)];
s2 = [s2 ;zeros(max_length-length(s2), 1)];

mix_signal = s1 + s2;




