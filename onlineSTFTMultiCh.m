clear;

resampFreq = 16e3;

[srcSig(:,:,1), sampFreq] = audioread('drums.wav'); % signal x channel x source (source image)
[srcSig(:,:,2), sampFreq] = audioread('piano.wav'); % signal x channel x source (source image)
srcSigResample(:,:,1) = resample(srcSig(:,:,1), resampFreq, sampFreq, 100); % resampling for reducing computational cost
srcSigResample(:,:,2) = resample(srcSig(:,:,2), resampFreq, sampFreq, 100); % resampling for reducing computational cost

% Mix source images of each channel to produce observed mixture signal
mixSig(1, :) = srcSigResample(:,1,1) + srcSigResample(:,1,2);
mixSig(2, :) = srcSigResample(:,2,1) + srcSigResample(:,2,2);
if abs(max(max(mixSig))) > 1 % check clipping
    error('Cliping detected while mixing.\n');
end

[nCh, nSamples] = size(mixSig);

nFft = 2048;
nOverlap = 2;
nDelay = nOverlap - 1;
nHop = nFft / nOverlap;
winFunc = hann(nFft, "periodic").';

nFrame = floor(nSamples / nHop) - 1;

sigDelay = zeros(nHop*nDelay, nCh).';
sigOut = zeros(nSamples, nCh).';

processedSig = zeros(nFft, nCh).';
processedSigDelay = zeros(nHop*nDelay, nCh).';

sigPartial = zeros(nFft, nCh);

scalingtmp = zeros(nFft*2, 1).';
scalingtmp(1:nFft) = winFunc.^2;
scalingtmp(nHop+1:nHop+nFft) = scalingtmp(nHop+1:nHop+nFft) + winFunc.^2;

scaling = 1./scalingtmp(nHop+1:nHop+nHop);



for iFrame=1:nFrame
    %% Extract Signal
    startIndex = (iFrame - 1) * nHop + 1;
    endIndex = startIndex + nHop - 1;    
    sig = mixSig(:, startIndex:endIndex);
   

    %% make FFT Frame
    sigPartial = [sigDelay sig];
    
    %% FFT
    spectra = zeros(nCh, nFft);
    for iCh=1:nCh
        spectra(iCh,:) = fft(sigPartial(iCh, :) .* winFunc);
    end
    spectraHalf = spectra(:, 1:nFft/2+1);
    
    %% processing
    magnitude = abs(spectraHalf);
    phase = angle(spectraHalf);


    processedSigSpectraHalf = magnitude .* exp(1j*phase);
    processedSigSpectra = [processedSigSpectraHalf flip(conj(processedSigSpectraHalf(:, 2:nFft/2)), 2)];
    
    %% IFFT
    for iCh=1:nCh
        processedSig(iCh,:) = real(ifft(processedSigSpectra(iCh,:))) .* winFunc;
    end
    
    %% Overlap Add
    sigOut(:, startIndex:endIndex) = (processedSig(:, 1:nHop) + processedSigDelay) .* scaling;


    %% Delay
    processedSigDelay = processedSig(:, nHop+1:nFft);
    sigDelay = sig;

end

plot([zeros(1, nHop) mixSig(1,:)]);
hold on;
plot(sigOut(1,:));
hold off;
legend(["original", "processed"]);
figure();
Expect = [zeros(nCh, nHop) mixSig];
Expect = Expect(1, 1:endIndex);
Result = sigOut(1, 1:endIndex);
plot(Expect - Result);