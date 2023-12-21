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

nFft = 512;
nOverlap = 2;
nDelay = nOverlap - 1;
nHop = nFft / nOverlap;
% winFunc = hann(nFft, "periodic").';
w = hann(nFft, "periodic").';
winFunc=sqrt(w*2.0*nHop/nFft);


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
for iCh=1:nCh
    for iBin=1:nBins
        Ukf(iCh,:,:,iBin) = eye(nCh) * 1e-3;
    end
end
UkfZ1= Ukf;

powerSpectra = zeros(nCh, nBins);
YkfEst = zeros(nCh, nBins);
forgetCoefs = 0.9;


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
    %% processing
    for iIter=1:1
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
%         r(iFrame,:) = sqrt(sum(powerSpectra, 2));
        r = sqrt(sum(powerSpectra, 2)); % (nSource x 1)
        %% Calculate weighted Covariance ( sphercical laplace distribution )
        weight = 1 ./ (2 * r) * 0.01;
%         weight = (1-forgetCoefs) * (1 ./ (r * 2));
        
        for iBin=1:nBins
            for iCh=1:nCh
                Ukf(iCh, :, :, iBin) = forgetCoefs * squeeze(Ukf(iCh, :, :, iBin)) + weight(iCh) * (Xkf(:, iBin) * Xkf(:, iBin)');
            end
        end

        %% Update demix matrix
        Wfnew = zeros(nCh);
        for iBin=1:nBins
            Wf = W(:, :, iBin);
            D = eye(nCh, nCh);
            for iCh=1:nCh
                WU = Wf * squeeze(Ukf(iCh, :, :, iBin)) + eye(nCh)*1e-6;
                WUinv = inv(WU);
                wkf = WUinv(:,iCh);

%                 d = sqrt(wkf' * squeeze(Ukf(iCh, :, :, iBin)) * wkf);
%                 wkf = wkf / d;
                Wfnew(:, iCh) = wkf;

                D(iCh,:) = wkf';
            end
            %% Scaling
            D = 
            W(:,:,iBin) = Wfnew;
%             A=inv(D);
% 
%             if abs(A(1, 1))>=abs(A(2, 1))
%                 a1=A(1, 1);
%             else
%                 a1=A(2, 1);
%             end
% 
%             if abs(A(2, 2))>=abs(A(1, 2))
%                 a2=A(2, 2);
%             else
%                 a2=A(1, 2);
%             end
% 
%             Wfnew=eye(nCh, nCh);
%             Wfnew(1, :)=a1*D(1,:)';
%             Wfnew(2, :)=a2*D(2,:)';
% 
%             W(:,:,iBin) = D;

        end
%         v = zeros(nCh, 1);  % source stering
%         for iBin=1:nBins
%             for iCh=1:nCh
%                 for iiCh=1:nCh
%                     d = W(:, iiCh, iBin)' * squeeze(Ukf(iCh,:,:, iBin)) * W(:, iiCh, iBin);
%                     if iCh~=iiCh
%                         u = W(:, iCh, iBin)' * squeeze(Ukf(iCh,:,:, iBin)) * W(:, iCh, iBin);
%                         v(iiCh, 1) = u / d;
%                     else
%                         v(iiCh, 1) = 1 - 1/sqrt(d);
%                     end
%                 end
%                 W(:,:, iBin) = W(:,:, iBin) - v * W(:, iCh, iBin)';
%             end
%         end
    end
%     cost(iFrame) = local_calcCostFunction(powerSpectra, W, nBins, 1)

    %% Bakc projection
%     eA = zeros(nCh);
%     for iBin=1:nBins
%         A = inv(W(:, :, iBin));
%         for iCh=1:nCh
%             eA(iCh,iCh) = A(1,iCh);
%         end
%         W(:,:,iBin) = eA * W(:,:, iBin);
%     end
%     D = zeros(nCh, nCh);
%     for iBin=1:nBins
%         for iCh=1:nCh
%             D(iCh,:) = W(:,iCh,iBin)';
%         end
%         A = inv(D) .* eye(nCh);
%         W(:,:,iBin) = D * W(:,:,iBin);
%     end
%         
    %% Separation
    Ykf = complex(zeros(nCh, nBins));
    for iBin=1:nBins
        Ykf(:, iBin) = W(:,:, iBin) * Xkf(:, iBin);
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