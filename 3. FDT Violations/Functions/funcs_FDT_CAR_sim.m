function [Csim,Asim,Rsim] = funcs_FDT_CAR_sim(tsaux,forceaux,noiseaux,NPARCELLS,Tup,Temp)
    %==========================================================================
    % FUNCTION: funcs_FDT_CAR_sim
    %--------------------------------------------------------------------------
    % Calculates the C, Asymmetry (A), and Response (R) functions
    % for a single simulation run. These metrics are used to test the 
    % Fluctuation-Dissipation Theorem (FDT).
    %
    % DEFINITIONS (based on Cugliandolo 1994):
    %   1. Correlation C(t,s): Autocorrelation of the state x.
    %      C(t,s) = x(t) * x(s)
    %
    %   2. Asymmetry A(t,s): Measures lack of detailed balance.
    %      A(t,s) = F(t)*x(s) - F(s)*x(t)
    %
    %   3. Response R(t,s): Linear response function using the FDT relation
    %      with respect to the noise term.
    %      R(t,s) = (1 / 2*Temp) * (x(t) * noise(s))
    %
    % INPUTS:
    %   ts          - [NPARCELLS x Tmax] Simulated time series (x)
    %   force       - [NPARCELLS x Tmax] Deterministic force term (-wC*x)
    %   noise       - [NPARCELLS x Tmax] Stochastic noise term (η)
    %   NPARCELLS   - Number of brain regions (nodes)
    %   Tmax        - Maximum time points (length of simulation)
    %   Temp        - System "Temperature", defined as sigma^2 / 2
    %
    % OUTPUTS:
    %   Csim        - [NPARCELLS x Tmax x Tmax] Correlation tensor
    %   Asim        - [NPARCELLS x Tmax x Tmax] Asymmetry tensor
    %   Rsim        - [NPARCELLS x Tmax x Tmax] Response tensor
    %
    % METHOD:
    %   Uses semi-vectorized operations over the time history (s < t) to 
    %   compute the outer products efficiently or nested 'ss' loops.
    %   Returns tensors with lower-triangular structure in time (causality).
    %
    % AUTHOR:
    %   Tomas Berjaga Buisan
    %==========================================================================

    % --- Memory Allocation ---
    % Pre-allocate 3D tensors: [Nodes x Time(t) x Time(s)]
    Csimaux = zeros(NPARCELLS,Tmax,Tmax);
    Asimaux = zeros(NPARCELLS,Tmax,Tmax);
    Rsimaux = zeros(NPARCELLS,Tmax,Tmax);

    % First option: nested loops
    for tt = 1:Tup
       for ss = 1:tt  % must be ss <= tt (causality)
           Csimaux(:,tt,ss) = tsaux(:,tt).*tsaux(:,ss);
           Asimaux(:,tt,ss) = forceaux(:,tt).*tsaux(:,ss) - forceaux(:,ss).*tsaux(:,tt);
           Rsimaux(:,tt,ss) = tsaux(:,tt).*noiseaux(:,ss);
       end
    end

    % Second option: matricial operations
    % for tt = 1:Tup
    %   tsaux_tt = tsaux(:, tt);
    %   forceaux_tt = forceaux(:, tt);
    %   Csimaux(:, tt, 1:tt) = bsxfun(@times, tsaux_tt, tsaux(:, 1:tt));
    %   Asimaux(:, tt, 1:tt) = bsxfun(@times, forceaux_tt, tsaux(:, 1:tt)) - bsxfun(@times, forceaux(:, 1:tt), tsaux_tt);
    %   nRsimaux(:, tt, 1:tt) = bsxfun(@times, tsaux_tt, noiseaux(:, 1:tt));
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --- Normalization ---
    % Normalize Response function by Temperature (1/sigma^2)
    Rsimaux = (0.5 / Temp) * Rsimaux;
end
