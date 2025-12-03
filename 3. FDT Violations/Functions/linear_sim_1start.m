function [ts, force, noise] = linear_sim_1start(dsig, dt, Tmax, TR, NPARCELLS, wC, x0)
    %==========================================================================
    % FUNCTION: linear_sim_1start
    %--------------------------------------------------------------------------
    % Simulates time series using a linear Langevin model, starting from
    % an initial thermalized state `x0`, using Euler–Maruyama integration.
    %
    % The system follows:
    %     dx/dt = -wC * x + η(t)
    % where:
    %     - wC is the effective connectivity (Jacobian),
    %     - η(t) is Gaussian white noise scaled by dsig.
    %
    % This function outputs the simulated time series (`ts`), along with the
    % corresponding force and noise at each sampled time point.
    %
    % INPUTS:
    %   dsig        - Noise amplitude scaled by sqrt(dt)
    %   dt          - Time step for simulation
    %   Tmax        - Number of time points to output (at TR resolution)
    %   TR          - Sampling interval for output (e.g., matching empirical data)
    %   NPARCELLS   - Number of brain regions (nodes)
    %   wC          - Effective connectivity matrix (G.C)
    %   x0          - Initial state vector (e.g., from thermalization)
    %
    % OUTPUTS:
    %   ts          - [NPARCELLS x Tmax] simulated time series
    %   force       - [NPARCELLS x Tmax] deterministic force (–wC*x)
    %   noise       - [NPARCELLS x Tmax] stochastic input (η), scaled by sqrt(dt)
    %
    % METHOD:
    %   The simulation runs at time resolution `dt`, but returns data sampled
    %   every `TR` seconds, aligned with BOLD or LFP sampling rates.
    %
    % AUTHOR:
    %   Tomas Berjaga Buisan
    %==========================================================================
    % Initialize state from provided x0
    x = x0;
    
    % Allocate memory for output
    xs     = zeros(Tmax, NPARCELLS);
    force  = zeros(Tmax, NPARCELLS);
    noise  = zeros(Tmax, NPARCELLS);

    % Total simulation time
    Tsim = (Tmax - 1) * TR;

    % Simulation counter for TR-aligned samples
    sample_idx = 0;
    % Integration loop
    for t = 0:dt:Tsim
        % Generate stochastic input (Gaussian noise)
        eta = dsig * randn(NPARCELLS, 1);
    
        % Compute deterministic force: –wC * x
        f = -wC * x;
    
        % Euler–Maruyama update
        x = x + dt * f + eta;
    
        % Sample system state every TR seconds
        if abs(mod(t, TR)) < dt / 100
            sample_idx = sample_idx + 1;
            xs(sample_idx, :)    = x(:)';
            force(sample_idx, :) = f(:)';
            noise(sample_idx, :) = eta(:)';
        end
    end

% Format outputs: [NPARCELLS x Tmax]
ts     = xs';
force  = force';
noise  = noise' / sqrt(dt);  % Normalize noise to match η(t) ~ √dt * ξ

end