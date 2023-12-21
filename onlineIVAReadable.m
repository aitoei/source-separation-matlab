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

nFft = 4096;
nOverlap = 2;
nDelay = nOverlap - 1;
nHop = nFft / nOverlap;
% winFunc = hamming(nFft, "periodic").';
w=hann(nFft, 'periodic');
winFunc=sqrt(w*2.0*nHop/nFft).';

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
nBins = nFft/2+1;

%% 
Wf = complex(eye(nCh));
W  = repmat(Wf, 1, 1, nBins); % (nCh x nCh x nBins)

est = complex(zeros(nBins, nCh));


Ukf  = zeros(nCh, nCh, nCh, nBins);
for iBin=1:nBins
    for iCh=1:nCh
        Ukf(iCh,:,:, iBin) = complex(eye(nCh))*1e-3;
    end
end

powerSpectra = zeros(nCh, nBins);
YkfEst = zeros(nCh, nBins);
forgetCoefs = 0.999;

isFirst = true;
for iFrame=1:nFrame
    %% Extract Signal
    startIndex = (iFrame - 1) * nHop + 1;
    endIndex = startIndex + nFft - 1;
    sig = mixSig(:, startIndex:endIndex);
   

    %% make FFT Frame
    sigPartial = sig;
    
    %% FFT
    spectra = zeros(nCh, nFft);
    for iCh=1:nCh
        spectra(iCh,:) = fft(sigPartial(iCh, :) .* winFunc);
    end
    spectraHalf = spectra(:, 1:nBins);
    
    Xkf = spectraHalf;
    Ykf = zeros(size(Xkf));
    %% processing
    for iIter=1:2
        %% Multiply demix matrix for each bin
        for iBin=1:(nBins)
            YkfEst(:, iBin) = W(:,:, iBin) * Xkf(:, iBin);
        end
        %% Calculate power spectrogram
        for iCh=1:nCh
            for iBin=1:nBins
                powerSpectra(iCh, iBin) = real(YkfEst(iCh, iBin)).^2 + imag(YkfEst(iCh, iBin)).^2;
            end
        end
        %% Calculate source model (Gauss) 全周波数が共通分散
        r = sum(powerSpectra, 2); % (nSource x 1)
        %% Calculate weighted Covariance ( sphercical laplace distribution )
%         weight = 1 ./ (2 * r);
        weight = (1 ./ (2 * sqrt(r)));
        for iBin=1:nBins
            %% Update covariance matrix
            for iCh=1:nCh
                Ucur = (weight(iCh) * Xkf(:, iBin) * Xkf(:, iBin)');
                Upre = squeeze(Ukf(iCh, :, :, iBin));
                U = Upre * forgetCoefs +  (1 - forgetCoefs) * Ucur;
                Ukf(iCh, :, :, iBin) = U;
            end
            %% Update demix matrix
            Wf = W(:, :, iBin)';
            Wnew = zeros(nCh, nCh);
            for iCh=1:nCh
                U = squeeze(Ukf(iCh, :, :, iBin));
                WU = Wf * U  + eye(nCh)*1e-6;
                WUinv = inv(WU);
                wi = WUinv(:, iCh);
                wi = wi / sqrt(wi' * U * wi);
                Wnew(:, iCh) = wi;
            end

            %% Projection back
            A = pinv(Wnew);
            
            scaleW = zeros(nCh, nCh);
            for iCh=1:nCh
                scaleW(iCh, iCh) = A(1, iCh);
            end
            
            %% Separation
            Ykf(:, iBin) = Wnew * Xkf(:, iBin);
            Ykf(1, iBin) = Ykf(1, iBin) * A(1,1);
            Ykf(2, iBin) = Ykf(2, iBin) * A(1,2);
            
            W(:, :, iBin) = Wnew;
        end
    end
    cost(iFrame) = local_calcCostFunction(powerSpectra, W, 1, nBins);
    processedSigSpectraHalf = Ykf;
    processedSigSpectra = [processedSigSpectraHalf flip(conj(processedSigSpectraHalf(:, 2:nFft/2)), 2)];
    
    %% IFFT
    for iCh=1:nCh
        processedSig(iCh,:) = real(ifft(processedSigSpectra(iCh,:))) .* winFunc;
    end
    
    %% Overlap Add
    sigOut(:, startIndex:endIndex) = sigOut(:, startIndex:endIndex) + (processedSig);
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

%% Local function for calculating cost function value in IVA
function [ cost ] = local_calcCostFunction(P, W, I, J)
logDetAbsW = zeros(I,1);
for i = 1:I
    logDetAbsW(i,1) = log(max(abs(det(W(:,:,i))), eps));
end
cost = sum(sqrt(sum(P, 3)), 'all')/J - 2*sum(logDetAbsW, 'all');
end