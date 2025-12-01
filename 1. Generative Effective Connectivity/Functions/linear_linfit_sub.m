function [Ceffsub,FCemp,FCsim] = linear_linfit_sub(tsdata,NPARCELLS,TR,Tinf,sigma,maxC, varargin)
% LINEAR_LINFIT_SUB  Fit an effective connectivity matrix C using FC(0) and lagged covariance
%                    for a single subject under a multivariate OU model.
%
%   [Ceff, FCemp, FCsim] = LINEAR_LINFIT_SUB(tsdata, nParcels, TR, Tinf, sigma, maxC, ...)
%
% Inputs
%   tsdata    : empirical time series.
%   nParcels  : number of parcels (integer).
%   TR        : sampling interval (seconds).
%   Tinf      : integer number of samples to trim before FC/COV.
%   sigma     : (scalar) noise variance for the OU model simulator.
%   maxC      : (scalar) normalization applied to C.
%
% Name-Value pairs (optional)
%   'Tau'        : lag in samples for lagged covariance target (default = 1).
%                  It can be modified depending on the data.
%   'EpsFC'      : learning rate for FC(0) term (default = 1e-4).
%                  MODIFY DEPENDING ON THE FITTING LEVELS
%   'EpsTau'     : learning rate for lagged covariance term (default = 1e-4).
%                  MODIFY DEPENDING ON THE FITTING LEVELS
%   'ErrorTol'   : stop when (old - new)/new < ErrorTol (default = 1e-5).
%                  Increase as possible...
%   'MaxIter'    : maximum number of iterations (default = 5000).
%   'UseAbsCov'  : use abs(cov) (true) to stabilize the variance normalization (default = true).
%
% Outputs
%   Ceff   : [nParcels x nParcels] fitted effective connectivity matrix.
%   FCemp  : [nParcels x nParcels] empirical FC(0) (Pearson correlation).
%   FCsim  : [nParcels x nParcels] simulated FC(0) from final model.
%
% Notes
%
%   • The update keeps C non-negative and re-enforces diagonal dominance to
%     satisfy a Hurwitz-like stability criterion at each step.

% Author  : Tomas Berjaga Buisan

% Copyright
%   © 2025 Your Lab. MIT License.
% ========================================================================

% ──────────────────────────────────────────────────────────────────────────
% Parse & validate inputs
% ──────────────────────────────────────────────────────────────────────────
p = inputParser;
p.addParameter('Tau',       1,       @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('EpsFC',     4e-4,    @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('EpsTau',    1e-4,    @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('ErrorTol',  1e-5,    @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MaxIter',   5000,    @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Seed',      [],      @(x)isempty(x)||(isnumeric(x)&&isscalar(x)));
p.addParameter('UseAbsCov', true,    @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});

Tau      = p.Results.Tau;
epsFC    = p.Results.EpsFC;
epsTau   = p.Results.EpsTau;
errTol   = p.Results.ErrorTol;
maxIter  = p.Results.MaxIter;
seed     = p.Results.Seed;
useAbs   = p.Results.UseAbsCov;

% ──────────────────────────────────────────────────────────────────────────
% Empirical targets: FC(0) and normalized lagged covariance COV_tau
% ──────────────────────────────────────────────────────────────────────────

indexN = 1:NPARCELLS;
N = length(indexN);
ts2  = tsdata(indexN,Tinf:end-Tinf);   % trim edges

% Empirical zero-lag FC (Pearson)
FCemp = corrcoef(ts2');   
% Empirical covariance for variance normalization
COVemp = cov(ts2');

if useAbs
    COVemp = abs(COVemp);
    FCemp = abs(FCemp);
end

% Lagged covariance at lag = Tau samples (normalized)
COVtauEmp = zeros(N, N);
tst = ts2';

for i = 1:N
    for j = 1:N
        sigratio(i,j) = 1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
        [clag, lags] = xcov(tst(:,i),tst(:,j),Tau);
        indx = find(lags == Tau);
        COVtauemp(i,j) = clag(indx)/size(tst,1);
     end
end

COVtauemp = COVtauemp.*sigratio;

% ──────────────────────────────────────────────────────────────────────────
% Initialize connectivity C (nonnegative, symmetric, diagonally dominant)
% ──────────────────────────────────────────────────────────────────────────

C = abs(randn(N));         % random nonnegative weights
C = (C + C.')/2;                  % symmetrize
diagonal  = sum(abs(C),2) - abs(diag(C));
C = diag(diagonal+eps) + C.*~eye(N);  % enforce diagonal dominance (Hurwitz-like)
%C = C + 0.1 * randn(NPARCELLS); % Add random asymmetry (OPTIONAL)
C = C/max(max(C))*maxC;        % normalize to maxC

% ──────────────────────────────────────────────────────────────────────────
% Optimize C by matching FC(0) and normalized COV(Tau)
% ──────────────────────────────────────────────────────────────────────────

Cnew = C;
olderror = inf;

for iter = 1:maxIter
    % --- Simulate OU / linear model ---
    [FCsim, COVsim, A] = linear_int(Cnew, sigma);

    % Compute normalized lagged covariance from model
    COVtausim = expm((Tau*TR)*A)*COVsim;
    COVtausim = COVtausim(1:N,1:N);
    for i = 1:N
        for j = 1:N
            sigratiosim(i,j) = 1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim = COVtausim.*sigratiosim;


    % Stopping criterion (checked every 10 iters)
    if mod(iter,10) < 0.1
        errornow = mean(mean((FCemp - FCsim).^2)) + mean(mean((COVtauemp - COVtausim).^2));
        if  (olderror - errornow)/errornow < errTol
            % disp('---> (olderror - errornow)/errornow < error_tol')
            break;
        end
        if  olderror < errornow
            % disp('---> olderror < errornow')
            break;
        end
        olderror = errornow;
    end

    %%% Learning - Gradient-like update
    for i = 1:N  
        for j = 1:N
            if (Cnew(i,j) > 0 || j == N-i+1)
                Cnew(i,j) = Cnew(i,j) + epsFC*(FCemp(i,j) - FCsim(i,j)) ...
                                      + epsTau*(COVtauemp(i,j) - COVtausim(i,j));
                if Cnew(i,j) < 0
                    Cnew(i,j) = 0;
                end
            end
        end
    end

    % Re-enforce diagonal dominance (Hurwitz-like) and renormalize
    diagonal = sum(abs(Cnew), 2) - abs(diag(Cnew));
    Cnew = diag(diagonal+eps) + Cnew.*~eye(N);  % enforce diagonal dominance (Hurwitz-like)
    Cnew = Cnew/max(max(Cnew))*maxC;
end
Ceffsub = Cnew;
end
