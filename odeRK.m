function sol = odeRK(ode, x0fcn, T, M,params)
   
    delta = 0.0001;  % Tolerance for convergence
    t = linspace(0, T, M+1);  % Time vector (ensure it's a row vector)
    dt = T / M;  % Time step size
    dt2 = dt / 2; 
    dt6 = dt / 6;
    n = numel(x0fcn([]));  % Number of state variables
    x = zeros(n, M+1);  % Solution array
    x(:, 1) = x0fcn([])';  % Initial conditions
    cache = zeros(n, 4);  % Cache for intermediate steps
    
    test = -1;
    while(test < 0)
        ox = x;  % Store the previous solution
        for i = 1:M
            % Runge-Kutta 4th-order method
            cache(:, 1) = ode(x(:, i),params);  % Evaluate the ODE at the current state
            cache(:, 2) = ode(x(:, i) + dt2 * cache(:, 1),params);  % Evaluate at intermediate step 2
            cache(:, 3) = ode(x(:, i) + dt2 * cache(:, 2),params);  % Evaluate at intermediate step 3
            cache(:, 4) = ode(x(:, i) + dt * cache(:, 3),params);  % Evaluate at intermediate step 4
            
            % Update the state using the Runge-Kutta method
            x(:, i + 1) = x(:, i) + dt6 * cache * [1 2 2 1]'; 
        end
        
        % Check for convergence: minimum error between previous and current solution
        test = min(delta * sum(abs(x), 2) - sum(abs(ox - x), 2));
    end
    sol.t = t';   % Time points (ensure it's a column vector)
    sol.y = x;    % Solution matrix (state variables at each time point)
    sol.solver = 'odeRK';  % Solver name (you can change it to your liking)
    sol.problem = 'Initial value problem';  % Problem description (can be modified)
    sol.stats.nfevals = M;  % Number of function evaluations (for tracking)
    sol.solver_type = 'Runge-Kutta 4th order';  % Type of solver used (you can modify it)
    
    % Create an interpolant for the solution using PCHIP (piecewise cubic Hermite interpolation)
    sol.interpolate = @(t_interp) pchip(sol.t, sol.y', t_interp);  % PCHIP interpolation (transposed solution matrix)
    sol.deval = @(t_eval) deval_compat(sol, t_eval);
end

function y_eval = deval_compat(sol, t_eval)
       y_eval = pchip(sol.t, sol.y, t_eval);  % No transpose of sol.y here
end
