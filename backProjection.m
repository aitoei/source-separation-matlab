clear;

[s1, fs1] = audioread("BASIC5000_0269.wav");
[s2, fs2] = audioread("TRAVEL1000_0282.wav");

assert(fs1 == fs2);
targetFs = 16e3;
[p, q] = rat(targetFs/fs1);
s1 = resample(s1, p, q);
s2 = resample(s2, p, q);

max_length = max(length(s1), length(s2));


s1 = [s1 ;zeros(max_length-length(s1), 1)];
s2 = [s2 ;zeros(max_length-length(s2), 1)];
original_signal  = [s1 s2].';
estimated_signal = original_signal .* 0.7;

mix_matrix = [ 0.9 0.5;
               0.3 0.7;
             ];


ovserved_signal = mix_matrix * original_signal;
ovserved_signal = ovserved_signal.';

referenceMic = 1;


[nSamples, nCh] = size(ovserved_signal);
nFft = 2048;
nOverlap = 2;
nHop = nFft / nOverlap;
winFunc = hann(nFft);

nFrame = floor(nSamples / nHop) - 1;
sigDelay = zeros(nHop, nCh);
sigOut = zeros(nSamples, nCh);
processedSig = zeros(nFft, nCh);
processedSigDelay = zeros(nHop, nCh);

for iFrame=1:nFrame
    %% Extract Signal
    startIndex = (iFrame - 1) * nHop + 1;
    endIndex = startIndex + nHop - 1;    
    ovserved_sig = ovserved_signal(startIndex:endIndex);
    sigPartial = [sigDelay; ovserved_sig];
    
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
    sigDelay = ovserved_sig;
    processedSigDelay = processedSig(nHop+1:end); 
end








