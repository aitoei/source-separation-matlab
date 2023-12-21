clear;

% Mix source images of each channel to produce observed mixture signal
[mixSig(1, :), fs] = audioread("mix_online_0.wav");
[mixSig(2, :), fs] = audioread("mix_online_1.wav");
[mixSig(3, :), fs] = audioread("mix_online_2.wav");
% if abs(max(max(mixSig))) > 1 % check clipping
%     error('Cliping detected while mixing.\n');
% end

[nCh, nSamples] = size(mixSig);

nFft = 4096;
nOverlap = 2;
nDelay = nOverlap - 1;
nHop = nFft / nOverlap;
winFunc = hamming(nFft, "periodic").';


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
        Ukf(iCh,:,:, iBin) = complex(eye(nCh))*1e-9;
    end
end

powerSpectra = zeros(nCh, nBins);
YkfEst = zeros(nCh, nBins);
forgetCoefs = 0.99;


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
    spectraHalf = spectra(:, 1:nBins);
    
    Xkf = spectraHalf;
    Ykf = zeros(size(Xkf));
    
    %% processing
    for iIter=1:5
        YkfEst = zeros(1, nBins);
        powerSpectra = zeros(1, nBins);
        for iCh=1:nCh
            for iBin=1:(nBins)
                YkfEst(iBin) = W(iCh,:, iBin) * Xkf(:, iBin);
            end

            for iBin=1:nBins
                powerSpectra(iBin) = real(YkfEst(iBin)).^2 + imag(YkfEst(iBin)).^2;
            end
            %% Calculate source model
            r = sqrt(sum(powerSpectra)); % (nCh x 1)
            %% Calculate contrast function
            weight = 1 ./ (2 * r); % (nCh x 1)
            %% Calculate covariance

            for iBin=1:nBins
                Ucur = weight * Xkf(:, iBin) * Xkf(:, iBin)';
                Upre = squeeze(Ukf(iCh, :, :, iBin));
                U = Upre * forgetCoefs +  (1 - forgetCoefs) * Ucur;
                Ukf(iCh, :, :, iBin) = U;
            end
            
            Wnew = zeros(nCh, nBins);
            E = eye(nCh);
            for iBin=1:nBins
                %% Update demix matrix
                Wf = W(:, :, iBin);
                U = squeeze(Ukf(iCh,:,:,iBin));
                WU = Wf  * U;
                WUinv = inv(WU);
                w = conj(WUinv(:, iCh));
                denom = sqrt(w' * U * w);
%                 w = conj(WU\E(:,iCh));
%                 denom = sqrt(w * U * w');
%                 w = w / sqrt(denom);
                Wnew(:,iBin) = w / denom;
            end
            
            W(iCh,:,:) = Wnew;
        end
    end

    %% restore scaling (projection back method)
    scaleW = zeros(nCh, nCh);
    for iBin=1:nBins
         A = inv(W(:,:, iBin));
         for iCh=1:nCh
             scaleW(iCh,iCh) = A(1,iCh);
         end
         W(:,:,iBin) = scaleW * W(:,:,iBin);
         Ykf(:, iBin) = W(:, :, iBin) * Xkf(:, iBin);
    end


    processedSigSpectraHalf = Ykf;
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
