function [FC, CV, A] = linear_int(gC, sigma)
    %% ========================================================================
    % FUNCTION: linear_int
    % -------------------------------------------------------------------------
    % Computes the steady-state covariance matrix (CV) and functional 
    % connectivity (FC) from a linear Langevin model with Jacobian matrix A = -gC.
    %
    % INPUTS:
    %   gC     - [NxN] Effective connectivity matrix
    %   sigma  - Scalar noise standard deviation
    %
    % OUTPUTS:
    %   FC     - Functional connectivity matrix (Pearson correlation)
    %   CV     - Covariance matrix obtained from Lyapunov equation
    %   A      - System Jacobian
    %
    % Author  : Tomas Berjaga Buisan
    % ========================================================================

    N = size(gC, 1);          % Number of brain regions

    % Construct Jacobian matrix: A = - gC
    A = - gC;

    % Input noise covariance (white, uncorrelated)
    Qn = (sigma^2) * eye(N);  

    % Solve the continuous Lyapunov equation: A*CV + CV*A' = -Qn
    CV = sylvester(A, A', -Qn);

    % Ensure symmetry and positivity
    CV = abs(CV + CV') / 2;

    % Convert to functional connectivity (Pearson correlations)
    FC = corrcov(CV);

end
