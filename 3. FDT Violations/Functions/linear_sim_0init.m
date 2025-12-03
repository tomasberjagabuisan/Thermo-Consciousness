function x0 = linear_sim_0init(dsig,dt,NPARCELLS,wC)
    %==========================================================================
    % FUNCTION: linear_sim_0init
    %--------------------------------------------------------------------------
    % Performs a thermal initialization (burn-in) for the linear Langevin 
    % model to stabilize the system dynamics before recording data.
    %
    % The system follows:
    %     dx/dt = -wC * x + η(t)
    %
    % This function runs the simulation for a fixed number of steps (2000)
    % and discards the trajectory, returning only the final state vector.
    % NOTE: Consider that if the system is complicated to converge to the
    % analyitical covariance. More burn-in will be needed.
    %
    % INPUTS:
    %   dsig        - Noise amplitude scaled by sqrt(dt)
    %   dt          - Time step for simulation
    %   NPARCELLS   - Number of brain regions (nodes)
    %   wC          - Effective connectivity matrix (Jacobian)
    %
    % OUTPUTS:
    %   x0          - [NPARCELLS x 1] Final state vector after burn-in
    %
    % METHOD:
    %   Initialize x with scaled random noise, then evolve using 
    %   Euler-Maruyama integration for 2000 time steps to reach a 
    %   stochastic steady state (thermalization).
    %
    % AUTHOR:
    %   Tomas Berjaga Buisan
    %==========================================================================
    
    % --- Configuration ---
    n_burnin = 5000;  % Number of steps to discard

    % --- Initialization ---
    % Initialize x scaled by the variance of noise (thermal initialization)
    % to ensure the system starts in the correct magnitude range.
    x = (dsig / sqrt(dt)) * randn(NPARCELLS, 1);

    % --- Burn-in Loop (Discard Transients) ---
    for step = 0:n_burnin
        % Generate stochastic input (Gaussian noise)
        eta = dsig * randn(NPARCELLS, 1);
        
        % Compute deterministic force: -wC * x
        f = -wC * x;
        
        % Euler-Maruyama update
        x = x + dt * f + eta;
    end
    % --- Final Output ---
    x0 = x;
end
