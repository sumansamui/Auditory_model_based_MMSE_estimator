% ----------------------------------------------------------------------
%        Auditory-model based Bayesian estimator for Speech enhancement
%            
%   This function implements an audio denoising algorithm based o paramaterized 
%   Bayesian estimator. The parameters The parameters of the estimator are chosen 
%   based on the human auditory characteristics such as the cochlea’s dynamic
%   compressive non-linearity and loudness perception theory. It uses an unbiased 
%   MMSE-Based Noise Power Estimation with Low Complexity and Low Tracking Delay
%   
%   This function takes noisy audio file as input and generates an enhanced audio file.  
%
%   Usage:  Auditory_MMSE_STSA_unbiased_ne(noisyFile,enhancedFile,alpha,beta)
%
%   Copyright (c) 2015 by Suman Samui
%
% ----------------------------------------------------------------------

function Auditory_MMSE_STSA_unbiased_ne(noisyFile,enhancedFile,alpha, beta)

[y, sr] = audioread(noisyFile);

[noisePowMat_p] = noisePowProposed(y,sr);

%% STFT parameters

frLen   = 32e-3*sr;  % frame size
fShift  = frLen/2;   % fShift size
nFrames = floor(length(y)/fShift)-1; % number of frames
anWin  = hanning(frLen,'periodic'); %analysis window
NFFT = 512;


%PERC=50; % window overlap in percent of frame size
len1=floor(fShift);
len2=frLen-len1;
%% Noise power calculation calculations - assuming that the first 6 frames is
% noise/silence 



NF = size(noisePowMat_p,2);

FB = size(noisePowMat_p,1);

q=1;
for z = 1:NFFT
    if z<=FB
       noisePowMat(z,:) =  noisePowMat_p(z,:);
        
    else
       noisePowMat(z,:) = noisePowMat_p(FB-q,:);
       q = q + 1;  
    end  
    
end

%% Processing
%--- allocate memory and initialize various variables
k = 1;
x_old=zeros(len1,1);
alpha1 = 0.9;
zeta_min=10^(-25/10);
eta=0.15; 
mu=0.98;


%alpha = 0.5;

%beta = 0.33;

%===============================  Start Processing =======================================================


for n=1:nFrames


    seg = anWin.*y(k:k+frLen-1);

    %--- Take fourier transform of  frame

    spec=fft(seg,NFFT);
    spec_ph = angle(spec);
    sig=abs(spec); % compute the magnitude
    sig2=sig.^2;
    
     
    % ----------------- estimate/update noise psd --------------
    noise_mu2 = noisePowMat(:,n);
    
    noise_mu=sqrt(noise_mu2);  % magnitude spectrum
    % ---------------------------------------------------------
    
    
    gammak = min(sig2./noise_mu2,40);  % posteriori SNR
    

    if n==1
        zeta = alpha1 + (1-alpha1)*max(gammak-1,0);
    else
        zeta = alpha1*Yk_prev./noise_mu2 + (1-alpha1)*max(gammak-1,0);     % a priori SNR
        zeta = max(zeta_min,zeta);  % limit zeta to -25 dB
    end
    
    
    log_sigma_k = gammak.* zeta./ (1+ zeta)- log(1+ zeta);    
    vad_decision= sum(log_sigma_k)/frLen;    
    if (vad_decision< eta) 
        % noise only frame found
        noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
    end
    % ===end of vad===
    
    A = zeta./(1+zeta);
        
    nu = A.*gammak ;
    
    C = sqrt(nu)./gammak ;
    
    N = gamma((beta/2)-alpha + 1).*confhyperg(alpha-(beta/2),1,-nu,100);
    
    D = gamma(-alpha + 1).*confhyperg(alpha,1,-nu,100);
    
    R = (N./D);
    
    R = R.^(1./beta);
    
    hw = C.*R;
        
    %hw = 1;

    sig=sig.*hw;
    Yk_prev=sig.^2;

    xi_w = ifft(hw.*spec, NFFT);
    xi_w = real(xi_w);


    % --- Overlap and add ---------------
    %
    xfinal(k:k+ len2-1)= x_old + xi_w(1:len1);
    x_old= xi_w(len1+ 1: frLen);

    k=k+len2;
end
%========================================================================================



audiowrite(enhancedFile,xfinal,sr);