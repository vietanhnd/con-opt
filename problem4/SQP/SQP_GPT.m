function [x_optimal, f_optimal, exit_flag, num_iterations, lambda, gradient, hessian] = SQP_GPT(objective, initial_guess, lower_bounds, upper_bounds, inequality_lower_bounds, inequality_upper_bounds, equality_constraint, inequality_constraint, method, hessian_func, max_iterations, tolerance)

    % Set default values
    if nargin < 12
        tolerance = 1e-3;
    end
    if nargin < 11
        max_iterations = 100;
    end
    if nargin < 10
        hessian_func = @dampedBFGS;
    end

    % Check feasibility of initial guess
    if any(initial_guess < lower_bounds) || any(initial_guess > upper_bounds)
        fprintf('Initial guess is not feasible. Adjusting...\n');
        initial_guess = (lower_bounds + upper_bounds) / 2;
    end

    % Initialization
    k = 0;
    converged = false;

    % Initial constraint values and gradients
    [equality_values, equality_gradient] = equality_constraint(initial_guess);
    [inequality_values, inequality_gradient] = inequality_constraint(initial_guess);
    inequality_values = [inequality_values - inequality_lower_bounds; inequality_upper_bounds - inequality_values];
    inequality_gradient = [inequality_gradient; -inequality_gradient];
    num_equality_constraints = length(equality_values);
    num_inequality_constraints = length(inequality_values);

    % Store information
    x_optimal = zeros(length(initial_guess), max_iterations);
    df = zeros(length(initial_guess), max_iterations);
    f = zeros(1, max_iterations);
    lambda_equality = zeros(num_equality_constraints, max_iterations);
    lambda_inequality = zeros(num_inequality_constraints, max_iterations);
    lambda_lower_bounds = zeros(length(initial_guess), max_iterations);
    lambda_upper_bounds = zeros(length(initial_guess), max_iterations);
    x_optimal(:, 1) = initial_guess;
    [f(1), df(:, 1)] = objective(initial_guess);

    % Initial Hessian
    H = zeros(length(initial_guess), length(initial_guess), max_iterations); 
    H(:, :, 1) = 10 * eye(length(initial_guess));

    % Main loop
    while k < max_iterations && ~converged
        k = k + 1;

        % Solve QP subproblem
        [p_k, ~, ~, ~, l] = quadprog(H(:, :, k), df(:, k), -inequality_gradient, inequality_values, equality_gradient', -equality_values, lower_bounds - x_optimal(:, k), upper_bounds - x_optimal(:, k));

        % Lagrange multiplier directions
        delta_y = l.eqlin - lambda_equality(:, k);
        delta_z = l.ineqlin - lambda_inequality(:, k);
        delta_a = l.lower - lambda_lower_bounds(:, k);
        delta_b = l.upper - lambda_upper_bounds(:, k);

        % Line search or trust region
        if strcmp(method, 'line-search')
            alpha = lineSearch(x_optimal(:, k), objective, equality_constraint, inequality_constraint, inequality_lower_bounds, inequality_upper_bounds, p_k, lambda_equality(:, k), lambda_inequality(:, k));
        elseif strcmp(method, 'trust-region')
            alpha = trustRegion(x_optimal(:, k), objective, equality_constraint, inequality_constraint, inequality_lower_bounds, inequality_upper_bounds, lambda_equality(:, k), lambda_inequality(:, k), p_k);
        else
            error('Invalid method. Please choose either line-search or trust-region');
        end

        % Update Lagrange multipliers
        lambda_equality(:, k+1) = lambda_equality(:, k) + alpha * delta_y;
        lambda_inequality(:, k+1) = lambda_inequality(:, k) + alpha * delta_z;
        lambda_lower_bounds(:, k+1) = lambda_lower_bounds(:, k) + alpha * delta_a;
        lambda_upper_bounds(:, k+1) = lambda_upper_bounds(:, k) + alpha * delta_b;

        % Compute updates
        D_x_Lagrange_1 = df(:, k) - equality_gradient * lambda_equality(:, k) - inequality_gradient' * lambda_inequality(:, k) - lambda_lower_bounds(:, k) + lambda_upper_bounds(:, k);
        x_optimal(:, k+1) = x_optimal(:, k) + alpha * p_k;
        [f(k+1), df(:, k+1)] = objective(x_optimal(:, k+1));
        [equality_values, equality_gradient] = equality_constraint(x_optimal(:, k+1));
        [inequality_values, inequality_gradient] = inequality_constraint(x_optimal(:, k+1));
        inequality_values = [inequality_values - inequality_lower_bounds; inequality_upper_bounds - inequality_values];
        inequality_gradient = [inequality_gradient; -inequality_gradient];
        D_x_Lagrange_2 = df(:, k+1) - equality_gradient * lambda_equality(:, k+1) - inequality_gradient' * lambda_inequality(:, k+1) - lambda_lower_bounds(:, k+1) + lambda_upper_bounds(:, k+1);
        p = x_optimal(:, k+1) - x_optimal(:, k);
        q = D_x_Lagrange_2 - D_x_Lagrange_1;

        % Update Hessian
        if nargin < 10
            H(:, :, k+1) = dampedBFGS(H(:, :, k), p, q);
        else
            H(:, :, k+1) = hessian_func(x_optimal(:, k+1), lambda_equality(:, k+1), lambda_inequality(:, k+1));
        end

        % Check convergence
        converged = norm(D_x_Lagrange_2, inf) < tolerance && ...
                    norm(p, inf) < tolerance && ...
                    norm(q, inf) < tolerance && ...
                    norm(delta_y, inf) < tolerance && ...
                    norm(delta_z, inf) < tolerance && ...
                    norm(delta_a, inf) < tolerance && ...
                    norm(delta_b, inf) < tolerance;
    end

    % Output
    f_optimal = f(k);
    lambda.eq = lambda_equality(:, k+1);
    lambda.ineq = lambda_inequality(:, k+1);
    lambda.lower = lambda_lower_bounds(:, k+1);
    lambda.upper = lambda_upper_bounds(:, k+1);
    gradient = df(:, k+1);
    hessian = H(:, :, k+1);
    num_iterations = k;

    % Exit flag
    if ~converged
        exit_flag = 0;
        fprintf('Maximum number of iterations reached without convergence.\n');
    else
        exit_flag = 1;
        fprintf('Successfully converged.\n');
    end
end
