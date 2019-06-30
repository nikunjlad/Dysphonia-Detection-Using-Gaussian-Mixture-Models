%% Without I-Vector feature extraction code 

close all;
clear;
clc;

%% DYSPHONIA DATABASE

Pathdir = uigetdir(cd,'Load directory for pathology recordings');
if ~isequal(Pathdir,0)
    d = dir(Pathdir);   % Use subdirectory 'wav' to read WAVE files
    d = d(3:end);
    dpath = cell(length(d), 1);
    k = 0;  % Loop index
    for i=1:length(d)
        [~,name,extn]=fileparts(d(i).name);
        if isequal(extn,'.wav')     % Select file only if valid WAVE file
            k = k + 1;
            % Read full address of WAVE file
            dpath{k} = fullfile(Pathdir, d(i).name);
        end
    end
    dpath = dpath(1:k, 1);
else
    msgbox('Select proper directory','Warning','warn')
end

%% NORMAL DATABASE

Pathdir = uigetdir(cd,'Load directory for pathology recordings');
if ~isequal(Pathdir,0)
    d = dir(Pathdir);   % Use subdirectory 'wav' to read WAVE files
    d = d(3:end);
    ppath = cell(length(d), 1);
    k = 0;  % Loop index
    for i=1:length(d)
        [~,name,extn]=fileparts(d(i).name);
        if isequal(extn,'.wav')     % Select file only if valid WAVE file
            k = k + 1;
            % Read full address of WAVE file
            ppath{k} = fullfile(Pathdir, d(i).name);
        end
    end
    ppath = ppath(1:k, 1);
else
    msgbox('Select proper directory','Warning','warn')
end

%% CONSTANTS

fs = 44000;
win = 0.020;
step = 0.010;

%% Feature Extraction

% Creating source and target variables
PF = cell(2, length(dpath));

tic;
%% Pathology recordings
for i = 1: length(dpath)
    s = audioread(dpath{i},[1 23000]);
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
    [jitter_local,jitter_local_absolute,jitter_RAP,jitter_PPQ5,shimmer_rel,...
        shimmer_local_absolute,shimmer_APQ3,shimmer_APQ5,shimmer_dB]=JitterShimmer(s,fs);
    
    combo = [EE E Z R C F mfcc lpcs lpccs rcs lars lsfs sscs formants ...
        jitter_local jitter_local_absolute jitter_RAP jitter_PPQ5 shimmer_rel...
        shimmer_local_absolute shimmer_APQ3 shimmer_APQ5];  % combo_dys = 102x51
    [combo_dys, ~, ~] = featureNormalize(combo);
    PF{1, i} = combo_dys;  % 2x50 cell 2 = no of classes (Dysphonia,Normal); 50 = number of .wav files/speakers
end

%% Normal recordings
for i = 1: length(ppath)
    s = audioread(ppath{i},[1 23000]);
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
    [jitter_local,jitter_local_absolute,jitter_RAP,jitter_PPQ5,shimmer_rel,...
        shimmer_local_absolute,shimmer_APQ3,shimmer_APQ5,shimmer_dB]=JitterShimmer(s,fs);
    
    combo = [EE E Z R C F mfcc lpcs lpccs rcs lars lsfs sscs formants ...
        jitter_local jitter_local_absolute jitter_RAP jitter_PPQ5 shimmer_rel...
        shimmer_local_absolute shimmer_APQ3 shimmer_APQ5];  % combo_dys = 102x51
    [combo_norm, ~, ~] = featureNormalize(combo);
    PF{2, i} = combo_norm;  % 2x50 cell 2 = no of classes (Dysphonia,Normal); 50 = number of .wav files/speakers
end
toc;

save PF.mat PF;

[classes,samples] = size(PF);

for i = 1:classes
    for j = 1:samples
        PF{i,j} = mean(PF{i,j});
    end
end

y = ones(samples * classes,1);
y(samples + 1:end,:) = 0;
PF = PF';
X = cell2mat(PF(:));
[X,~] = featureNormalize(X);
[Xtr,ytr,Xcv,ycv,Xtst,ytst,X,y] = randomize(X,y);

save dt Xtr ytr Xcv ycv Xtst ytst X y;