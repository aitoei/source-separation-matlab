clear;


filename = "BASIC5000_0269.wav";
[signal, fs] = audioread(filename);
signal = signal.';

[nCh, nSamples] = size(signal);

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

rng(10, "twister");
wn = dsp.ColoredNoise("Color","white", "BoundedOutput",true, "SamplesPerFrame", nHop, "RandomStream", "Global stream");
noisePower = 0.1;
noiseMag = zeros(nFft/2+1, nCh).';

Input = zeros(size(signal));

for iFrame=1:nFrame
    %% Extract Signal
    startIndex = (iFrame - 1) * nHop + 1;
    endIndex = startIndex + nHop - 1;    
    sig = signal(:, startIndex:endIndex);
    noise = (wn() * noisePower).';
    sig = sig + noise;
    Input(startIndex:endIndex) = sig;
   

    %% make FFT Frame
    sigPartial = [sigDelay sig];
    
    %% FFT
    sigPartial = sigPartial .* winFunc;    
    spectra = fft(sigPartial);
    spectraHalf = spectra(1:nFft/2+1);
    
    %% processing
    %% Noise Estimation
    if iFrame <= 2
        noiseMag = noiseMag + abs(spectraHalf)./2;
    end

    magnitude = abs(spectraHalf);
    phase = angle(spectraHalf);

    %% Spectral Subtraction
    for iBin=1:(nFft/2)+1
        magnitude(iBin) = magnitude(iBin) - noiseMag(iBin);
        if magnitude(iBin) < 0
            magnitude(iBin) = 0;
        end
    end


    processedSigSpectraHalf = magnitude .* exp(1j*phase);
    processedSigSpectra = [processedSigSpectraHalf flip(conj(processedSigSpectraHalf(2:nFft/2)))];
    
    %% IFFT
    processedSig = real(ifft(processedSigSpectra)) .* winFunc;
    
    %% Overlap Add
    sigOut(:, startIndex:endIndex) = (processedSig(1:nHop) + processedSigDelay) .* scaling;


    %% Delay
    processedSigDelay = processedSig(:, nHop+1:nFft);
    sigDelay = sig;

end

plot([zeros(1, nHop) Input(1,:)]);
hold on;
plot(sigOut(1,:));
hold off;
legend(["original", "processed"]);
