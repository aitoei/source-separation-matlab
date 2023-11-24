clear;


filename = "BASIC5000_0269.wav";
[signal, fs] = audioread(filename);

[nSamples, nCh] = size(signal);

nFft = 2048;
nOverlap = 4;
nDelay = nOverlap - 1;
nHop = nFft / nOverlap;
winFunc = hann(nFft);

nFrame = floor(nSamples / nHop) - 1;


sigDelay = zeros(nHop*nDelay, nCh);
sigDelayPtr = 0;
sigOut = zeros(nSamples, nCh);

processedSig = zeros(nFft, nCh);
processedSigDelay = zeros(nHop*nDelay, nCh);
processedSigDelayPtr = 0;

sigPartial = zeros(nFft, nCh);
getBuffer = zeros(nHop, nCh);

for iFrame=1:nFrame
    %% Extract Signal
    startIndex = (iFrame - 1) * nHop + 1;
    endIndex = startIndex + nHop - 1;    
    sig = signal(startIndex:endIndex);

    %% Get delay signal from Overlap Buffer
    sigPartial(1:nDelay*nHop) = sigDelay;
    sigPartial(nDelay*nHop+1:nFft) = sig;
    
    %% Overlap Buffer Update
    sigDelay(1:nHop*(nDelay-1)) = sigDelay(nHop+1:nHop*nDelay);
    sigDelay(nHop*(nDelay-1)+1:nHop*nDelay) = sig;
    
    %% FFT
    sigPartial = sigPartial .* winFunc;    
    spectra = fft(sigPartial);
    
    %% processing
    % do nothing
    
    %% IFFT
    processedSig = ifft(spectra) .* winFunc;
    
    %% Overlap Add
    sigOut(startIndex:endIndex) = processedSig(1:nHop) + processedSigDelay(1:nHop);

    %% Buffer Update
    processedSigDelay(1:nHop*(nDelay-1)) = processedSigDelay(nHop+1:nHop*nDelay);
    processedSigDelay(nHop*(nDelay-1)+1:nHop*nDelay) = 0;

    %% Overlap Add
    processedSigDelay = processedSigDelay + processedSig(nHop+1:nFft);

end
