function [noisePowMat] = noisePowProposed(noisy,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% propose SPP algorithm to estimate the spectral noise power
%%%% papers: "Unbiased MMSE-Based Noise Power Estimation with Low Complexity and Low Tracking Delay", IEEE TASL, 2012 
%%%% "Noise Power Estimation Based on the Probability of Speech Presence", Timo Gerkmann and Richard Hendriks, WASPAA 2011
%%%% Input :        noisy:  noisy signal
%%%%                   fs:  sampling frequency
%%%%                   
%%%% Output:  noisePowMat:  matrix with estimated noise power for each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Copyright (c) 2011, Timo Gerkmann
%%%%%%%%%%%%%%%%%%%%%% Author: Timo Gerkmann and Richard Hendriks
%%%%%%%%%%%%%%%%%%%%%% Universitaet Oldenburg
%%%%%%%%%%%%%%%%%%%%%% KTH Royal Institute of Technology
%%%%%%%%%%%%%%%%%%%%%% Delft university of Technology
%%%%%%%%%%%%%%%%%%%%%% Contact: timo.gerkmann@uni-oldenburg.de
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Universitaet Oldenburg, Delft university, KTH 
%	Royal Institute of Technology nor the names of its contributors may be 
% 	used to endorse or promote products derived from this software without 
% 	specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% Version v1.0 (October 2011)


%% some constants
frLen   = 32e-3*fs;  % frame size
fShift  = frLen/2;   % fShift size
nFrames = floor(length(noisy)/fShift)-1; % number of frames

anWin  = hanning(frLen,'periodic'); %analysis window

%% allocate some memory
noisePowMat = zeros(frLen/2+1,nFrames);


%% initialize
noisePow = init_noise_tracker_ideal_vad(noisy,frLen,frLen,fShift, anWin); % This function computes the initial noise PSD estimate. It is assumed that the first 5 time-frames are noise-only.
noisePowMat(:,1)=noisePow; 

PH1mean  = 0.5;
alphaPH1mean = 0.9;
alphaPSD = 0.8;

%% constants for a posteriori SPP
q          = 0.5; % a priori probability of speech presence:
priorFact  = q./(1-q);
xiOptDb    = 15; % optimal fixed a priori SNR for SPP estimation
xiOpt      = 10.^(xiOptDb./10);
logGLRFact = log(1./(1+xiOpt));
GLRexp     = xiOpt./(1+xiOpt);

for indFr = 1:nFrames
    indices       = (indFr-1)*fShift+1:(indFr-1)*fShift+frLen;
    noisy_frame   = anWin.*noisy(indices);
    noisyDftFrame = fft(noisy_frame,frLen);
    noisyDftFrame = noisyDftFrame(1:frLen/2+1);
	
    noisyPer = noisyDftFrame.*conj(noisyDftFrame);
    snrPost1 =  noisyPer./(noisePow);% a posteriori SNR based on old noise power estimate


    %% noise power estimation
	GLR     = priorFact .* exp(min(logGLRFact + GLRexp.*snrPost1,200));
	PH1     = GLR./(1+GLR); % a posteriori speech presence probability

	PH1mean  = alphaPH1mean * PH1mean + (1-alphaPH1mean) * PH1;
	stuckInd = PH1mean > 0.99;
	PH1(stuckInd) = min(PH1(stuckInd),0.99);
	estimate =  PH1 .* noisePow + (1-PH1) .* noisyPer ;
	noisePow = alphaPSD *noisePow+(1-alphaPSD)*estimate;
        
	noisePowMat(:,indFr) = noisePow;
end
return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function   noise_psd_init = init_noise_tracker_ideal_vad(noisy,fr_size,fft_size,hop,sq_hann_window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%This m-file computes an initial noise PSD estimate by means of a
%%%%Bartlett estimate.
%%%%Input parameters:   noisy:          noisy signal
%%%%                    fr_size:        frame size
%%%%                    fft_size:       fft size
%%%%                    hop:            hop size of frame
%%%%                    sq_hann_window: analysis window
%%%%Output parameters:  noise_psd_init: initial noise PSD estimate
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Author: Richard C. Hendriks, 15/4/2010
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
 
for I=1:5
    noisy_frame=sq_hann_window.*noisy((I-1)*hop+1:(I-1)*hop+fr_size);
    noisy_dft_frame_matrix(:,I)=fft(noisy_frame,fft_size);
end
noise_psd_init=mean(abs(noisy_dft_frame_matrix(1:fr_size/2+1,1:end)).^2,2);%%%compute the initialisation of the noise tracking algorithms.
return


% EOF
