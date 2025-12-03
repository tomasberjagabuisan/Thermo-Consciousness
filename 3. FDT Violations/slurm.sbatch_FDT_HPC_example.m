#!/bin/bash
#SBATCH --job-name=FDT   
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-1000
#SBATCH --output=FDT%A_%a.out
#SBATCH --error=FDT%A_%a.err

# Load MATLAB module
ml MATLAB

matlab -nojvm -nodisplay<<-EOF

%% ========================================================================
%  FDT SIMULATIONS
%  ------------------------------------------------------------------------
%  This script is intended to be run as an array job on a SLURM cluster.
%  Each array index "s" produces one independent simulation batch with
%  different random seeds, and saves partial results to disk. 
%
%  NOTE: You have to adjust the simulation time according to the relaxation
%  time of your matrix, please be sure to see first the ACF to set it
%  properly. Optionally for demanding jobs the GEC could be multiply by a 
%  constant the decay to go faster as an approximation.
%
%  NOTE 2: A faster analytical way is currently under development so please
%  contact the corresponding author beacuse maybe it is available at that
%  time!
% ========================================================================

%% 0. SLURM ARRAY INDEX & RANDOM SEED HANDLING
s = str2double(getenv('SLURM_ARRAY_TASK_ID'));

% Pause time (in seconds) to decorrelate RNG seeding across array tasks
pause_time = mod(s, 60);
pause(pause_time);

% Different seed every time (dependent on system time)
rng('shuffle');

%% 1. PATHS & BASIC CONFIGURATION
addpath('Functions');

% -------------------------------------------------------------------------
% DATASET SELECTION
%   DATASET = 1 --> LFP ANESTHESIA
%   DATASET = 2 --> EEG ANESTHESIA
%   DATASET = 3 --> EEG DoC
% -------------------------------------------------------------------------
DATASET     = 2;

% -------------------------------------------------------------------------
% CONDITION (used for DATASET 2, kept here for completeness)
%   1 --> WAKE EC
%   2 --> WAKE EO
%   3 --> XEN
%   4 --> PRO
%   5 --> KET
% -------------------------------------------------------------------------
CONDITION = 1;

% Frequency band label
Frequencies = '05_40';

% -------------------------------------------------------------------------
% MODEL / FITTING PARAMETERS
% -------------------------------------------------------------------------
SETSIGMA = 0.2;   % Std. dev. of noise for Ceff fitting 
                  % (Set lower if the matrix is big...)
maxC     = 0.2;   % Maximum value for the C matrix (default)
NSUBSIM  = 10000; % Number of simulations per subject (as big as necessary 
                  % to converge to more stable estimates). Used as in Monti
                  % et al., Phys Rev E, 2025


% Simulation initialization mode:
%   SIMINIT = 1 --> initializes every simulation independently
%   SIMINIT = 0 --> initializes only the first simulation (per subject),
%                    reuse the same z0 for subsequent simulations
SIMINIT    = 0;


% Optional preprocessing of simulated signals
DETRENDSIM = 0;   % 1: detrend after demeaning; 0: no detrend
ZSCORESIM  = 0;   % 1: z-score simulated TS; 0: raw TS

%% 2. LOAD EMPIRICAL DATA SETTINGS
% -------------------------------------------------------------------------
% Channels and sampling parameters
% -------------------------------------------------------------------------
myChannels = 1:60;
NPARCELLS  = numel(myChannels);

if strcmp(Frequencies, '05_4')
    Fs = 8;
else
    Fs = 80;
end

TR   = 1 / Fs;
Tmax = 1000; % NOTE: REALLY IMPORTANT TO SET IT PROPERLY PLEASE LOOK AT THE 
            % COMMENT ABOVED
disp(['Tmax = ', num2str(Tmax)])

Tinf = 1;
Tsup = Tmax - Tinf;
Tup = Tsup - Tinf + 1;

%% 3. LOAD EFFECTIVE CONNECTIVITY (Ceff)
% -------------------------------------------------------------------------
fileCeff = sprintf('Fittings/%s/CeffWAKEEC_%s.mat', Frequencies, Frequencies);
Ceffsub  = load(fileCeff);
NSUB = 10;   % number of subjects for this condition


%% 4. LINEAR MODEL – SETUP
% -------------------------------------------------------------------------
disp('-- 3. LINEAR Model --');

% -------------------------------------------------------------------------
% MODEL PARAMETERS
% -------------------------------------------------------------------------
dt = 0.1*TR/2; % Time step
sig = SETSIGMA; 
dsig = sqrt(dt)*sig; % noice variance
Temp = sig^2/2; % "Temperature"

% -------------------------------------------------------------------------
% Allocate subject-level containers
% -------------------------------------------------------------------------
FCsim_sub = zeros(NSUB,NPARCELLS,NPARCELLS);

Csub = zeros(NSUB,NPARCELLS,Tup,Tup);
Asub = zeros(NSUB,NPARCELLS,Tup,Tup);
Rsub = zeros(NSUB,NPARCELLS,Tup,Tup);
xRsub = zeros(NSUB,NPARCELLS,Tup,Tup);

dtCsub = zeros(NSUB,NPARCELLS,Tup,Tup);
dsCsub = zeros(NSUB,NPARCELLS,Tup,Tup);

iVFDTsub = zeros(NSUB,NPARCELLS,Tup,Tup);
xiVFDTsub = zeros(NSUB,NPARCELLS,Tup,Tup);

% -------------------------------------------------------------------------
% Check for partial results / define output file
% -------------------------------------------------------------------------
substart     = 1;
file_partial = sprintf('./Results/%s/WAKEEC_%s_%d.mat', Frequencies, Frequencies, s);

%% 5. MAIN LOOP OVER SUBJECTS
% -------------------------------------------------------------------------
for sub = substart:NSUB

    fprintf('+: Subject = %d of %d\n', sub, NSUB);

    % wC matrix (G.Cij) obtained by fitting FC and COVtau
    wC = cell2mat(squeeze(Ceffsub.CeffWAKEEC(sub)));
    
    % ---------------------------------------------------------------------
    % Initialize simulation-level containers
    % ---------------------------------------------------------------------
    tssim = zeros(NPARCELLS,Tup);
    fosim = zeros(NPARCELLS,Tup);
    nosim = zeros(NPARCELLS,Tup);

    FCsim = zeros(NPARCELLS,NPARCELLS);
    Csim = zeros(NPARCELLS,Tup,Tup);
    Asim = zeros(NPARCELLS,Tup,Tup);
    Rsim = zeros(NPARCELLS,Tup,Tup);

    % 5.1 LOOP OVER SIMULATIONS
    %  --------------------------------------------------------------------
    for sim = 1:NSUBSIM

        % Initialization discarding first 5000t time steps:
        % Either:
        %   - initialize every simulation (SIMINIT == 1)
        %   - or only the first one, then reuse x0 (SIMINIT == 0 && sim == 1)
        if SIMINIT == 1 || (SIMINIT == 0 && sim == 1)
            x0 = linear_sim_0init(dsig, dt, NPARCELLS, wC);
        end

        % Actual simulation
        [ts, force, noise] = linear_sim_1start(dsig,dt,Tmax,TR,NPARCELLS,wC,x0);

        % -----------------------------------------------------------------
        % Optional preprocessing of simulated signals (The mean can also be
        % removed and maybe it facilitates some calculations!)
        % -----------------------------------------------------------------
        if DETRENDSIM == 1
            for np = 1:NPARCELLS
                ts(np,:)    = detrend(demean(ts(np,:)));
                force(np,:) = detrend(demean(force(np,:)));
                noise(np,:) = detrend(demean(noise(np,:)));
            end
        end

        if ZSCORESIM == 1
            for np = 1:NPARCELLS
                ts(np,:)    = zscore(ts(np,:));
                force(np,:) = zscore(force(np,:));
                noise(np,:) = zscore(noise(np,:));
            end
        end

        % Restrict to analysis time window
        ts = ts(:,Tinf:Tsup);
        force = - force(:,Tinf:Tsup);
        noise = noise(:,Tinf:Tsup);
        
        % -----------------------------------------------------------------
        % FC simulated
        % -----------------------------------------------------------------
        FCsim = FCsim + corrcoef(ts');

        % -----------------------------------------------------------------
        % Compute C(t,s), A(t,s) and R(t,s) for this simulation
        % -----------------------------------------------------------------
        [Csimaux,Asimaux,Rsimaux] = funcs_FDT_CAR_sim(ts,force,noise,NPARCELLS,Tup,Temp);

        % Accumulate over simulations
        Csim = Csim + Csimaux;
        Asim = Asim + Asimaux;
        Rsim = Rsim + Rsimaux;
    end

    %% 5.2 AVERAGING OVER SIMULATIONS
    FCsim_sub(sub,:,:) = FCsim./NSUBSIM;
    clear FCsim

    Csub(sub,:,:,:) = Csim./NSUBSIM;
    Asub(sub,:,:,:) = Asim./NSUBSIM;
    Rsub(sub,:,:,:) = Rsim./NSUBSIM;

    clear Csim Asim Rsim

    %% 5.3 DIFFERENTIATION OF C IN t AND s
    dtCaux = zeros(NPARCELLS,Tup,Tup);
    dsCaux = zeros(NPARCELLS,Tup,Tup);

    dxaux = TR*(1:Tup);
    for np = 1:NPARCELLS
        Caux = squeeze(Csub(sub,np,:,:));

        for ti = 1:Tup
            dtCaux(np,:,ti) = derivative(dxaux,Caux(:,ti));
            dsCaux(np,ti,:) = derivative(dxaux,Caux(ti,:));

            % Uncommento to check differentiation: 
            % Ctest_t(np,:,ti) = cumtrapz(dxaux,squeeze(dtCaux(np,:,ti))) + Caux(1,ti);
            % Ctest_s(np,ti,:) = cumtrapz(dxaux,squeeze(dsCaux(np,ti,:))) + Caux(ti,1);
        end
    end

    dtCsub(sub,:,:,:) = dtCaux;
    dsCsub(sub,:,:,:) = dsCaux;

    % Uncommento to check differentiation: %
    % Ctest_t_sub(sub,:,:,:) = Ctest_t;
    % Ctest_s_sub(sub,:,:,:) = Ctest_s;
    % plot(squeeze(Csub(1,1,:,1)),'LineWidth',1.5), hold on,
    % plot(squeeze(Ctest_t_sub(1,1,:,1)),'LineWidth',1.5),
    % plot(squeeze(Ctest_s_sub(1,1,:,1)),'--','LineWidth',1.5), hold off

    clear Caux dtCaux dsCaux
    %% FDT METRICS
    % R(t,s) as in Cugliandolo 1994 eq. (2.9)
    xRsub = (1.0/2.0/Temp) * (dsCsub - dtCsub - Asub);

    % FDR as in Lippiello 1999 eq. (43)
    xXFDRsub = 0.5 * (1 - dtCsub ./ dsCsub) - 0.5 * Asub ./ dsCsub;

    % I remove NANs that may appear when dividing by zero-values of dsCsim
    xXFDRsub(isnan(xXFDRsub)) = 0;

    % OBS: also could be calculated as:
    % X = Temp * xRsub / dsCsub
    %   = Temp * Rsub / dsCsub
    XFDRsub = Temp * Rsub ./ dsCsub;

    % dVFDT Differential Violation of FDT
    % as defined in Cugliandolo 1997 eq. (1)
    % dVFDT = dC(t,s)/ds - T * R(t,s)
    dVFDTsub = dsCsub - Temp * Rsub;
    xdVFDTsub = dsCsub - Temp * xRsub;

    % iVFDT Integral violation of FDT
    % as defined in Cugliandolo 1997 eq. (2)
    iVFDTsub(sub,:,1,1) = 0;
    xiVFDTsub(sub,:,1,1) = 0;
    iVFDTsub(sub,:,Tup,Tup) = 0;
    xiVFDTsub(sub,:,Tup,Tup) = 0;

    for np = 1:NPARCELLS
        for tt = 2:Tup
            for ss = 1:tt-1
                tintaux = TR * (ss:tt);

                % R integral with original R
                Rintaux = squeeze(Rsub(sub,np,tt,(ss:tt)));
                intRaux = trapz(tintaux,Rintaux);
                iVFDTsub(sub,np,tt,ss) = Csub(sub,np,tt,tt) - Csub(sub,np,tt,ss) - Temp * intRaux;

                % R integral with xR
                Rintaux = squeeze(xRsub(sub,np,tt,(ss:tt)));
                intRaux = trapz(tintaux,Rintaux);
                xiVFDTsub(sub,np,tt,ss) = Csub(sub,np,tt,tt) - Csub(sub,np,tt,ss) - Temp * intRaux;
            end
        end
    end

    clear tintaux intRaux

    %% 5.8 SAVE PARTIAL PROGRESS
    sublast = sub;

    saveOK = 0;
    save(file_partial,'saveOK')
    save(file_partial,'sublast','-append')

    save(file_partial,'saveOK','Csub','Asub','Rsub','xRsub','dtCsub','dsCsub',...
                      'XFDRsub','xXFDRsub','dVFDTsub','xdVFDTsub','iVFDTsub','xiVFDTsub',...
                      'Temp','sublast', 'ts','noise','force','-v7.3') 

    saveOK = 1;
    save(file_partial,'saveOK','-append')
end

EOF
