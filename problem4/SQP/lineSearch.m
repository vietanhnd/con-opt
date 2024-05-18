function alpha = lineSearch(x_k, objective, h, g,g_l, g_u, p_k, y, z)
    %{
        Implements line search algorithm for SQP based on slide 25, Lec SQP.

        Uses Powell's l1 merit function: P(x, λ, μ) = f(x) + λ|h(x)| + μ|min {0, g(x)}|

        Inputs:
            x_k: Current iterate
            objective: Objective function handle (returns [f, df])
            h: Equality constraint function handle (returns [h, dh])
            g: Inequality constraint function handle (returns [g, dg])
                - Converted to standard form: G(x) = [g(x) - g_l; g_u - g(x)] etc
            g_l: Lower bound for inequality constraints
            g_u: Upper bound for inequality constraints
            p_k: Step direction
        Outputs:
            alpha: Optimal step size
    %}

    fprintf('------ Performing Line Search ---------\n');

    alpha = 1;
    stop = false;
    nEq = length(h(x_k));
   

    [f_k, df_k] = objective(x_k);
    [h_k, ~] = h(x_k);
    [g_k, ~] = g(x_k);
    G_k = [g_k - g_l; g_u - g_k];
    nIneq = length(G_k);
    
    % Penalty parameters
    % Penalty function: P(x, λ, μ) = f(x) + λ|h(x)| + μ|min {0, g(x)}|, where g(x) = [g(x) - g_l; g_u - g(x)] 
    % Penalty parameters: λ >= |y|, μ >= |z|, where y are Lagrange multipliers for equality constraints and z are Lagrange multipliers for inequality constraints

    lambda = ones(nEq, 1)*max(abs(y));
    mu = ones(nIneq ,1)*max(abs(z));

    % Merit function
    fprintf('|h_k|:\n');
    disp(abs(h_k));
    fprintf('|min {0, G_k}|:\n');
    disp(abs(min(0, G_k)));

    % c = φ(0) = f(xk) + λ|h(xk)| + μ|min {0, g(xk)|
    c = f_k + lambda'*abs(h_k) + mu'*abs(min(0, G_k));
    % b = φ′(0) = ∇f (xk)′∆xk − penalty
    b = df_k'*p_k - lambda'*abs(h_k) - mu'*abs(min(0, G_k));


    fprintf('b:\n');
    disp(b);
    fprintf('c:\n');
    disp(c);
   
    % Backtracking line search
    while ~stop

        fprintf('Alpha:\n');
        disp(alpha);

        x = x_k + alpha*p_k;

        [f_x, ~] = objective(x);
        [h_x, ~] = h(x);
        [g_x, ~] = g(x);
        G_x = [g_x - g_l; g_u - g_x];
        
        phi_alpha = f_x + lambda'*abs(h_x) + mu'*abs(min(0, G_x));


        if phi_alpha <= c + 0.1*alpha*b % Armijo condition
            fprintf('------ Armijo Condition Satisfied ---------\n');
            fprintf('New Alpha: %d\n', alpha);
            stop = true;
        else % Update alpha
            fprintf('update alpha\n');
            a = (phi_alpha-(c+b*alpha))/(alpha^2);
            fprintf('a:\n');
            disp(a);
            alpha_min = -b/(2*a);
            fprintf('alpha_min:\n');
            disp(alpha_min);
            alpha = min(0.9*alpha, max(0.1*alpha, alpha_min));
        end
    end
end