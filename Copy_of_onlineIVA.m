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
UkfZ1= zeros(nCh, nCh, nCh, nBins);

powerSpectra = zeros(nCh, nBins);
YkfEst = zeros(nCh, nBins);
forgetCoefs = 0.5;



nummics = 2;
numrefs = 0;

%% perform stft
addpath('stft2');
% fft size
fftsize=512;
stftshift=fftsize/2;

M=nummics;
R=numrefs;
N=M;


%% params
% forgetting factor for bss
BF_FORGET=0.999;
% the shape parameter of the source prior
GAMMA=0.2;
%
% used to keep stable
%
VAR_BIAS=0.01;
STABLE_EPS=1e-3;
BF_DIAGLOAD=1e-6;

%% space for bss
% data buffer
fsize=M+R;
Bufdata=zeros(nBins,fsize);
% current mic backup
Miccurrent=zeros(nBins,M);

% the weighted correlation matrices
C1=cell(nBins, 1);
C2=cell(nBins, 1);
for k=1:nBins
    C1{k}=STABLE_EPS*eye(fsize, fsize);
    C2{k}=STABLE_EPS*eye(fsize, fsize);
end

% demixing matrices
Demix=cell(nBins, 1);
for k=1:nBins
    Demix{k}=eye(fsize, fsize);
end

%% perform iteration
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

    
    %
    % calculate nonlinearity
    %
    phi1=0;
    phi2=0;
    
    for k=1:nBins
        x=Xkf(k,:).';
        y=Demix{k}*x;
        % output data
%         Bssout(nBins,:)=y(1:M).';
        
        phi1=phi1+abs(y(1))^2;
        phi2=phi2+abs(y(2))^2;
    end
    
    phi1 = (1-BF_FORGET)*(1 ./ (2 * sqrt(phi1)));
    phi2 = (1-BF_FORGET)*(1 ./ (2 * sqrt(phi2)));
%     phi1=(1-BF_FORGET)*(phi1+VAR_BIAS)^((GAMMA-2)/2);
%     phi2=(1-BF_FORGET)*(phi2+VAR_BIAS)^((GAMMA-2)/2);

    % update the demixing matrices
    for k=1:K
        %
        % accumulate the weighted correlation
        %
        x=Bufdata(nBins,:).';
        C1{k}=BF_FORGET*C1{k}+phi1*(x*x');
        C2{k}=BF_FORGET*C2{k}+phi2*(x*x');
        
        %
        % update demixing matrix
        %
        H1=inv(Demix{k}*C1{k}+BF_DIAGLOAD*eye(fsize, fsize));
        w1=H1(:, 1);
        
        H2=inv(Demix{k}*C2{k}+BF_DIAGLOAD*eye(fsize, fsize));
        w2=H2(:, 2);
        
        D=eye(fsize, fsize);
        D(1, :)=w1';
        D(2, :)=w2';
        
        %
        % solve the scaling ambiguity
        %
        A=inv(D);
        
        if abs(A(1, 1))>=abs(A(2, 1))
            a1=A(1, 1);
        else
            a1=A(2, 1);
        end
        
        if abs(A(2, 2))>=abs(A(1, 2))
            a2=A(2, 2);
        else
            a2=A(1, 2);
        end
        
        D=eye(fsize, fsize);
        D(1, :)=a1*w1';
        D(2, :)=a2*w2';
        
        Demix{k}=D;
    end
    
    for m=1:M
        Ytf{m}(:, iFrame)=Bssout(:, m);
    end
end

%% perform istft and output signal
dataout=zeros(dataLength(T, stftshift, fftsize ), N);
for n=1:N
    dataout(:, n)=istft(Ytf{n}, stftshift, false);
end