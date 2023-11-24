clear;


filename = "BASIC5000_0269.wav";
[signal, fs] = audioread(filename);

[nSamples, nCh] = size(signal);

nFft = 2048;
nOverlap = 2;
nHop = nFft / nOverlap;
winFunc = hamming(nFft, "periodic");

nFrame = floor(nSamples / nHop) - nOverlap;


sigDelay = zeros(nHop, nCh);
sigOut = zeros(nSamples, nCh);
processedSig = zeros(nFft, nCh);
processedSigDelay = zeros(nHop, nCh);

wSum = zeros(nSamples, nCh);
for iFrame=1:nFrame
    %% Extract Signal
    startIndex = (iFrame - 1) * nHop + 1;
    endIndex = startIndex + nFft - 1;    
    sig = signal(startIndex:endIndex);

    
    %% FFT
    sigPartial = sig .* winFunc;
    spectra = fft(sigPartial);
    
    %% processing
    % do nothing
    
    %% IFFT
    processedSig = (real(ifft(spectra)) .* winFunc);
    %% Overlap add
    sigOut(startIndex:endIndex) = sigOut(startIndex:endIndex) + processedSig;
    wSum(startIndex:endIndex) = wSum(startIndex:endIndex) + winFunc.^2;
end
pos = (wSum ~= 0);
sigOut(pos) = sigOut(pos) ./ wSum(pos);
plot(signal);
hold on;
plot(sigOut)
legend(["original", "fft"])