load('optimized.mat');
load('idx69.mat');
load('ubm.mat');

%% Feature Extraction
[Fname, Pname]= uigetfile('*.wav','Pick file for testing', cd);
s = audioread(fullfile(Pname,Fname), [1 23000]);

fs = 44000;            % Sampling frequency
win = 0.020;           % window length
step = 0.010;          % overlap length 

pf = cell(1,1);

tic;                   % Start stopwatch timer
% Pathology recordings
s = bsxfun(@minus,s,mean(s)) / max(abs(s));
EE= (Energy_Entropy_Block(s, win*fs, step*fs, 10))';
E = (ShortTimeEnergy(s, win*fs, step*fs));
Z = (zcr(s, win*fs, step*fs, fs));
R = (SpectralRollOff(s, win*fs, step*fs, 0.80, fs))';
C = (SpectralCentroid(s, win*fs, step*fs, fs));
F = (SpectralFlux(s, win*fs, step*fs, fs));
mfcc = (msf_mfcc(s,fs,'nfilt',40,'ncep',12));
lpcs = (msf_lpc(s,fs,'order',12));
lpccs = (msf_lpcc(s,fs,'order',12));
rcs = (msf_rc(s,fs,'order',12));
lars = (msf_lar(s,fs,'order',12));
lsfs = msf_lsf(s,fs,'order',12);
sscs = msf_ssc(s,fs,'nfilt',12);
formants = formant(s,fs);
[jitter_local,jitter_local_absolute,jitter_RAP,jitter_PPQ5,shimmer_rel,shimmer_local_absolute,...
    shimmer_APQ3,shimmer_APQ5,shimmer_dB]=JitterShimmer(s,fs);
combo_dys = [EE E Z R C F mfcc lpcs lpccs rcs lars lsfs sscs formants jitter_local...
    jitter_local_absolute jitter_RAP jitter_PPQ5 shimmer_rel...
    shimmer_local_absolute shimmer_APQ3 shimmer_APQ5];  % combo_dys = 102x51
[X_norm, ~, ~] = featureNormalize(combo_dys);
pf{1, 1} = X_norm;  % 2x50 cell 2 = no of classes (Dysphonia,Normal); 50 = number of .wav files/speakers
toc;

%% UBM Mapping to the incoming .wav file

nWorkers = 2;
nmix = 256;                             % In this case, we know the # of mixtures needed
final_niter = 10;                       % Number of iterations required
ds_factor = 1;                          % Down-Sampling Factor

las1 = cell2mat(pf(:));
[nspeaker, nchannel] = size(pf);         % # of classes, # of samples/class
lass1 = las1(:,idx);                     % new data <==> to ith feature set
m2c = mat2cell(lass1,ones(1,nchannel) .* 51,length(idx)); % convert lass to cell again
pf = reshape(m2c,[nspeaker,nchannel]);   % reshape the vector cell into matrix cell
% ubm = gmm_em(pf(:), nmix, final_niter, ds_factor, nWorkers);

stats1 = cell(1,size(pf,2));             % statistics required for I-Vector model

[N,F] = compute_bw_stats(pf{1,1}, ubm);
stats1{1,1} = [N; F];

tvDim = no_of_ivectors;           % 150
niter = 5;               
TT = train_tv_space(stats1(:), ubm, tvDim, niter, nWorkers);  % supervector

devIVs1 = zeros(tvDim, 1, size(pf,2));
devIVs1(:, 1, 1) = extract_ivector(stats1{1, 1}, ubm, TT);

save devIVsingle devIVs1;

%% Applying Cholesky whitening

[ivect,class,samples] = size(devIVs1);
tt = zeros(ivect,1);   
tt(:,1) = devIVs1(:,1);

alpha = 0.1;                               % Cholesky whitening co-efficient
W = zeros(size(tt,2));
W = W + cov(tt,1);

W = (1 - alpha) * W + alpha * eye(length(W));
vwccn = chol((W) ^ -1,'lower');
ivect = tt * vwccn;