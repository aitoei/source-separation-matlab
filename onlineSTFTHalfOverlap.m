clear;


filename = "BASIC5000_0269.wav";
[signal, fs] = audioread(filename);

[nSamples, nCh] = size(signal);

nFft = 2048;
nOverlap = 4;
nHop = nFft / nOverlap;
winFunc = hann(nFft);

nFrame = floor(nSamples / nHop) - 1;

nDelay = nOverlap - 1;
sigDelay = zeros(nHop*nDelay, nCh);
sigOut = zeros(nSamples, nCh);
processedSig = zeros(nFft, nCh);
processedSigDelay = zeros(nHop*nDelay, nCh);

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
    processedSig = ifft(spectra);
    
    %% Overlap add
    sigOut(startIndex:endIndex) = (processedSig(1:nHop) + processedSigDelay);
    
    %% Delay Update
    for iDelay=nDelay:1
        startIndex = (iDelay - 1) * nHop + 1;
        endIndex   = startIndex + nHop;
        sigDelay(startIndex+nHop:endIndex) = sigDelay(startIndex:endIndex);
    end
    sigDelay(startIndex:endIndex) = sig;

    sigDelay = sig;
    processedSigDelay = processedSig(nHop+1:end);    
end
