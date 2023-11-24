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
    for iCh=1:nCh
        for iSample=(0:nHop*nDelay-1)
            index = mod(iSample + sigDelayPtr, nHop*nDelay);
            sigPartial(iSample+1, iCh) = sigDelay(index+1, iCh);
        end
    end
    for iCh=1:nCh
        for iSample=(nHop*nDelay+1):nFft
            sigPartial(iSample) = sig(iSample-nHop*nDelay);
        end
    end
    
    %% Overlap Buffer Update
    for iCh=1:nCh
        for iSample=0:nHop-1
            index = mod(iSample + sigDelayPtr, nHop*nDelay);
            sigDelay(index+1) = sig(iSample+1);
        end
    end
    sigDelayPtr = mod(sigDelayPtr + nHop, nFft);

    
    %% FFT
    sigPartial = sigPartial .* winFunc;    
    spectra = fft(sigPartial);
    
    %% processing
    % do nothing
    
    %% IFFT
    processedSig = ifft(spectra) .* winFunc;

    %% Get delay processed signal from buffer
    for iSample=0:nHop-1
        index = mod(iSample + processedSigDelayPtr, nHop*nDelay);
        for iCh=1:nCh
            getBuffer(iSample+1) = processedSigDelay(index+1, iCh);
        end
    end
    processedSigDelayPtr = mod(processedSigDelayPtr + nHop, nFft);

    %% Overlap Add
    sigOut(startIndex:endIndex) = processedSig(1:nHop) + getBuffer;

    %% Update Processed buffer
    for iCh=1:nCh
        for iSample=(0:nHop*(nDelay-1)-1)
            index = mod(iSample + processedSigDelayPtr, nHop*nDelay);
            processedSigDelay(index+1, iCh) = processedSigDelay(index+1, iCh) + processedSig(nHop+iSample+1, iCh);
        end

        for iSample=nHop*(nDelay-1):nHop*nDelay-1
            index = mod(iSample + processedSigDelayPtr, nHop*nDelay);
            processedSigDelay(index+1, iCh) = processedSig(nHop+iSample+1, iCh);            
        end
    end
end
