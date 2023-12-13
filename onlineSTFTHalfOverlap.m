clear;


filename = "BASIC5000_0269.wav";
[signal, fs] = audioread(filename);

[nSamples, nCh] = size(signal);

nFft = 2048;
nOverlap = 2;
nHop = nFft / nOverlap;
winFunc = hann(nFft, "periodic");

nFrame = floor(nSamples / nHop) - 1;


sigDelay = zeros(nHop, nCh);
sigOut = zeros(nSamples, nCh);
processedSig = zeros(nFft, nCh);
processedSigDelay = zeros(nHop, nCh);

tmp = zeros(nFft*2, 1);
tmp(1:nFft) = winFunc.^2;
tmp(nHop+1:nHop+nFft) = tmp(nHop+1:nHop+nFft) + winFunc.^2;
tmp(nFft+1:nFft+nFft) = tmp(nFft+1:nFft+nFft) + winFunc.^2;
scaling = 1./ tmp(nHop+1:nHop+nFft);

for iFrame=1:nFrame
    %% Extract Signal
    startIndex = (iFrame - 1) * nHop + 1;
    endIndex = startIndex + nHop - 1;    
    sig = signal(startIndex:endIndex);
    sigPartial = [sigDelay; sig];
    
    %% FFT
    sigPartial = sigPartial .* winFunc;    
    spectra = fft(sigPartial);
    
    %% processing
    % do nothing
    
    %% IFFT
    processedSig = (real(ifft(spectra)) .* winFunc);
    
    %% Overlap add
    sigOut(startIndex:endIndex) = (processedSig(1:nHop) + processedSigDelay);
    
    %% Delay Update
    sigDelay = sig;
    processedSigDelay = processedSig(nHop+1:end);    
end
