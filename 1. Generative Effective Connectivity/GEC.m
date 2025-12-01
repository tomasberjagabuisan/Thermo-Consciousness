%% FITTING of EEG data by a MOU (Synthetic Data Example)
%
% ──────────────────────────────────────────────────────────────────────────
% WHAT THIS SCRIPT DOES
% ──────────────────────────────────────────────────────────────────────────
% This is a self-contained demonstration of Multivariate OU fitting
% pipeline without loading any external EEG files. It:
%   1) Generates synthetic multivariate time series X(t) 
%   2) Preprocesses the data as follows (modify if necessary):
%        • 0.5–40 Hz band filtering (low-pass @ 40 Hz, high-pass @ 0.5 Hz)
%        • downsampling (to ~80 Hz; TR = 1/80 s)
%        • detrending + z-scoring
%        • envelope extraction (upper envelope) - If desired
%   3) Fits an effective connectivity matrix Ceff by matching:
%        • empirical FC(0) (zero-lag Pearson correlation)
%        • normalized lagged covariance at lag Tau
%   4) Saves the averaged Ceff and diagnostics
%
% OUTPUTS
%   • Ceffsub  : fitted effective connectivity (averaged over NSIMFIT runs)
%   • FCempsub : empirical FC(0)
%   • FCsimsub : FC(0) from the final model

% Author  : Tomas Berjaga Buisan

% Copyright
%   © 2025 Your Lab. MIT License.
% ========================================================================

clear all; clc;
addpath Functions\

%%% SET DATASET  - If you have multiple DATASETS, you can label it
DATASET = 1;
% If you reuse this script inside a multi-condition repo layout, CONDITION
% can label “state”.
CONDITION = 2;

Frequencies = '05_40';

% ──────────────────────────────────────────────────────────────────────────
% Synthetic generation (Replace by Loading your Actual Time Series)
% ──────────────────────────────────────────────────────────────────────────

fs_native       = 400;     % native simulation rate (Hz) for OU data
dur_sec         = 400;     % total duration (s)
Tnative         = dur_sec * fs_native;
nParcels        = 60;      % channels/parcels

sparsity  = 0.10;          % fraction of non-zero off-diagonals
maxCtrue  = 0.2;           % cap for connectivity weights
tau_decay = 1.0;           % intrinsic decay time constant (s)
sigma_syn = 0.02;          % process noise for OU simulation


% Build a random, nonnegative symmetric Ctrue with given sparsity and cap.
Ctrue = zeros(nParcels);
mask  = triu(rand(nParcels) < sparsity, 1);
W     = max(0, randn(nParcels)); W = triu(W,1);
Ctrue(mask) = W(mask);
Ctrue = Ctrue + Ctrue.';                        % symmetrize

if max(Ctrue(:)) > 0, Ctrue = Ctrue ./ max(Ctrue(:)) * maxCtrue; end
diagDom = sum(abs(Ctrue),2) - abs(diag(Ctrue)); % enforce diagonal dominance
Ctrue   = diag(diagDom) + Ctrue.*~eye(nParcels);

% Stable OU drift A = -I/tau + C, then push eigenvalues left if needed.
Atrue = -eye(nParcels)/tau_decay + Ctrue;
eigA = eig(Atrue); maxReal = max(real(eigA)); alpha = 0.2;
if maxReal >= -alpha
    Atrue = Atrue - (maxReal + alpha) * eye(nParcels);
end

% Simulate OU via Euler–Maruyama: dX = A X dt + sqrt(2*sigma) dW
dt = 1/fs_native;
sqrtQ = sqrt(2*sigma_syn*dt);
burn = round(5*fs_native); % burn-in to reach stationarity
X = zeros(nParcels, Tnative);
x = zeros(nParcels,1);
for t = 1:(Tnative+burn)
    x = x + Atrue*x*dt + sqrtQ*randn(nParcels,1);
    if t > burn
        X(:, t-burn) = x;
    end
end

% ──────────────────────────────────────────────────────────────────────────
% Pre-processing (replaces for desired, here an example)
% ──────────────────────────────────────────────────────────────────────────
% Band of interest (mirrors your EEG preprocessing)
low_cutoff_frequency  = 0.5;  % high-pass (Hz)
high_cutoff_frequency = 40.0;  % low-pass  (Hz) 

DETRENDEMP   = 1;      % detrend empirical (here: synthetic) 
ZSCOREEMP    = 1;      % z-score empirical (here: synthetic) 
ENV          = 0;      % do envelope  empirical if needed (here: synthetic)

% NOTE: The envelope applies a non-linear transformation to the data, 
% so it may not be appropriate for fitting a multivariate OU process in some cases. 
% Keep this in mind when interpreting the results.

%% 1. Preprocess the data (0.5–40 Hz band → downsample → detrend/z-score)
% IMPORTANT: Downsampling is chosen conservatively as fs/(2*fc_high) to
% reduce aliasing risk after low-pass at 40 Hz.
downsampling_factor = max(1, round(fs_native/(2*high_cutoff_frequency)));

% Low-pass 40 Hz (order 3 to avoid ill-conditioned systems)
low_order = 3; % It can be changed!
[low_b, low_a] = butter(low_order, high_cutoff_frequency / (fs_native/2), 'low');

% High-pass 0.5 Hz
high_order = 3; % It can be changed!
[high_b, high_a] = butter(high_order, low_cutoff_frequency / (fs_native/2), 'high');

% Apply low-pass filter and then high-pass filter to each row (channel) of the data matrix
filtered_data_low  = zeros(size(X));
filtered_data_band = zeros(size(X));

for i = 1:size(X, 1)
    filtered_data_low(i, :) = filtfilt(low_b, low_a, X(i, :));
    filtered_data_band(i, :) = filtfilt(high_b, high_a, filtered_data_low(i, :));
    downsampled_data(i,:) = downsample(filtered_data_band(i, :), downsampling_factor);
end

tsdata = downsampled_data;

%%% Preprocessing: remove mean, detrend, and z-score.
%%% Z-scoring and detrending do not affect correlations,
%%% but help improve numerical stability and result accuracy.

EEGdata       = tsdata;
for np = 1:nParcels
    if DETRENDEMP == 1
       EEGdata(np,:) = EEGdata(np,:) - mean(EEGdata(np,:));
       EEGdata(np,:) = detrend(EEGdata(np,:));
    end
    if ZSCOREEMP == 1
            EEGdata(np,:) = zscore(EEGdata(np,:));
    end
    if ENV == 1
            [yupper,ylower] = envelope(EEGdata(np,:));
            EEGdata(np,:) = yupper;
    end 
end
tsdata = EEGdata;



%% 2: Fitting to obtain G.C (Effective Connectivity)

NSIMFIT  = 1; % Number of times the fitting is done (SET TO 1 FOR RUNNING FAST)
% NOTE: Is important to do a high time NSIMFIT as the inizialization is
% randomm. RECOMMENDED: 1000-10000
TR = 1/(2*high_cutoff_frequency); % Effective TR after downsampling
Tinf = 10; % point trimming
SETSIGMA_fit = 0.05;   % noise variance passed to linear_int during fitting
% NOTE: It can be adjusted depending on the problem and also made 
% heterogeneous if needed for studying focal problems
maxC         = 0.2;    % cap for C during iterative updates

disp('-- 2. Fitting to obtain Ceff = G.Cij --')

%%% Check if fitting was previously done

fileCeff = sprintf(['Ceff_COND_%d_NSIMFIT_%d_%s.mat'],...
                                 CONDITION,NSIMFIT,Frequencies);
if isfile(fileCeff)
    load(fileCeff)
    dolinfit = 0;
    disp('! Load Ceff from file')
else
    dolinfit = 1;
    sigma = SETSIGMA_fit;
end

%%% Fit of G.C using a Linear Model
if dolinfit == 1
    tsaux = tsdata;
    disp('! Do Fitting to obtain Ceff')
    Ceffsub = zeros(nParcels,nParcels);
    FCempsub = zeros(nParcels,nParcels);
    FCsimsub = zeros(nParcels,nParcels);

    for simfit = 1:NSIMFIT
        simfit;
        %disp(['+: Linear fit for sub. = ', num2str(sub), ' of ', num2str(NSUB)])
        [Ceffaux,FCempaux,FCsimaux] = linear_linfit_sub(tsaux,nParcels,TR,Tinf,sigma,maxC);
        Ceffsub = squeeze(Ceffsub(:,:)) + Ceffaux;
        FCempsub = squeeze(FCempsub(:,:)) + FCempaux;
        FCsimsub = squeeze(FCsimsub(:,:)) + FCsimaux;
    end

    clear tsaux Ceffaux Ceffgaux FCempaux FCsimaux
    Ceffsub = Ceffsub./NSIMFIT;
    FCempsub = FCempsub./NSIMFIT;
    FCsimsub = FCsimsub./NSIMFIT;
    Isubdiag = find(tril(ones(60),-1));
    cc=corrcoef(atanh(FCempsub(Isubdiag)),atanh(FCsimsub(Isubdiag)));
    corrFC = cc(2);
    %%% Save the Ceffsub to use later
    save(fileCeff,'Ceffsub','FCempsub','FCsimsub','corrFC')
end
