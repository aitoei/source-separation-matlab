clear;

filename = "BASIC5000_0269.wav";
[signal, fs] = audioread(filename);

nFft = 2048;
nOverlap = 2;
nHop = nFft / nOverlap;
winFunc = hann(nFft, "periodic");

iscola(winFunc, nHop, "ola")

s = stft(signal, fs, "Window", winFunc, "OverlapLength", nHop, "FFTLength", nFft);
[X,T] = istft(s,fs,'Window',winFunc,'OverlapLength',nHop,'FFTLength',nFft,'Method','ola','ConjugateSymmetric',true);

plot(signal);
hold on;
plot(X);