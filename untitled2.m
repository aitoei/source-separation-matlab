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
    for iIter=1:2
        for iCh=1:nCh
            y = zeros(1, nBins);
            for iBin=1:(nBins)
                y(iBin) = W(iCh,:, iBin) * Xkf(:, iBin);
            end
            powerSpectra = zeros(1, nBins);
            for iBin=1:nBins
                powerSpectra(iBin) = real(y(iBin)).^2 + imag(y(iBin)).^2;
            end
            %% Calculate source model
            r = sqrt(sum(powerSpectra)); % (1 x 1)
            %% Calculate contrast function
            weight = 1 ./ (2 * r); % (nCh x 1)
            %% 
            wnew = zeros(nCh, 1, nBins);
            for iBin=1:nBins
                %% Calculate covariance
                Ucur = (weight * Xkf(:, iBin)) * Xkf(:, iBin)';
                Upre = squeeze(Ukf(iCh, :, :, iBin));
                U = Upre * forgetCoefs +  (1 - forgetCoefs) * Ucur;
                Ukf(iCh, :, :, iBin) = U;
                %% Update demix matrix
                WU = Wf  * U;
                WUinv = inv(WU);
                w = WUinv(iCh, :);
                denom = w * U * w';
                w = w / sqrt(denom);
                wnew(:,1,iBin) = w';
            end
            W(iCh,:, :) = wnew;
            A = W(:,:,3);
        end
    end

    %% restore scaling (projection back method)
    eA = zeros(nCh, nCh);
    for iBin=1:nBins
         A = inv(W(:,:, iBin));
         for iCh=1:nCh
             eA(iCh,iCh) = A(1,iCh);
         end
         W(:,:,iBin) = eA * W(:,:,iBin);
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