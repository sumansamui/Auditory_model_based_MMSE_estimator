
addpath(genpath(pwd));

file = 'si1581';

filename.clean = [file '.wav']; 

sr = 16000; % experiment sampling rate

snr = 5;

filename.noise = 'babble_16k.wav';

filename.noisy = [file 'noisy_snr' num2str(snr) '_' filename.noise];
        
filename.enhanced = [file '_enhanced_proposed.wav'];
       
addnoise_asl(filename.clean, filename.noise, filename.noisy , snr);
        
[alpha, beta] = alpha_beta_selection(sr,512);

Auditory_MMSE_STSA_unbiased_ne(filename.noisy ,filename.enhanced, alpha, beta);

%MMSE_LSA(filename.noisy ,filename.enhanced);

pesq_noisy = pesq(filename.clean, filename.noisy);
        
pesq_enhanced = pesq(filename.clean, filename.enhanced);


%--------------------------------------------------------
enhanced = audioread(filename.enhanced);
        
clean = audioread(filename.clean);
        
noisy = audioread(filename.noisy);


stoi_noisy = stoi(clean,noisy,sr);

stoi_enhanced = stoi(clean,enhanced,sr);


%--------------------------------------------------------

