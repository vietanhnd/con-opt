%{
    Function performing SQP with line search for a general constrained optimisation problem.

    min f(x)
    s.t: 
        h(x) = 0
        g_l <= g(x) <= g_u
        x_l <= x <= x_u
    
    which has been converted to the form:

    min f(x)
    s.t:
        h(x) = 0
        g(x) - g_l >= 0
        g_u - g(x) >= 0
        x_l - x >= 0
        x - x_u >= 0

    Lagrange multipliers represented by:
    - y: Lagrange multipliers for equality constraints
    - z: Lagrange multipliers for inequality constraints
    - a: Lagrange multipliers for lower box constraints
    - b: Lagrange multipliers for upper box constraints


    Parameters:
        - f: function handle to objective function
        - x0: initial guess
        - x_l, x_u: Lower and upper bounds (box constraints)
        - g_l, g_u: Lower and upper bounds for inequality constraints
        - h: function handle to equality constraint
            - Returns [h, dh] where dh is the jacobian of h
        - g: function handle to inequality constraint
            - Returns [g, dg] where dg is the jacobian of g
        - method: 'line-search' or 'trust-region'   
            - Uses respective methods for determining step size.
        - hessianfunc (optional): User defined hessian of lagrange. (default uses damped BFGS approximation)
        - MaxIter (optional): Maximum number of iterations (default 100)
        - tol (optional): Tolerance for stopping criteria (default 1e-6)

    Returns:
        - x: array containing optimal solution at each iteration
        - fval: array containing function values at each iteration        
        - exitflag: exit flag
            - 1: converged
            - 0: maximum number of iterations reached/error
        - iterations: number of iterations
        - lambda: Struct containing lagrange multipliers
            - lambda.eq: Lagrange multipliers for equality constraints
            - lambda.ineq: Lagrange multipliers for inequality constraints
            - lambda.lower: Lagrange multipliers for lower box constraint
            - lambda.upper: Lagrange multipliers for upper box constraint
        - grad: array containing gradients at each iteration
        - hessian: array containing hessians at each iteration
%}


function [x, fval, exitflag, iterations, lambda, grad, H] = SQP(objective, x0, x_l, x_u, g_l, g_u, h, g, method, hessianfunc, MaxIter, tol)

    % Set default values
    if nargin < 12
        tol = 1e-6;
    end
    if nargin < 11
        MaxIter = 10;
    end
    if nargin < 10 % User has not specified hessian function
        hessianfunc = @dampedBFGS;
    end

    % Check if x0 is within bounds - if not place in the middle of the box
    if any(x0 < x_l) || any(x0 > x_u)
        fprintf('Initial guess is not feasible. Placing it in the middle of the box.\n');
        x0 = (x_l + x_u)/2;
    end


    n = length(x0);

    % Initialisation
    k = 0;
    converged = false;
    % Initial Constraint values and gradients
    [h_k, dh_k] = h(x0);
    % Writing inequality constraints in standard form
    [g_k, dg_k] = g(x0);
    G_k = [g_k - g_l; g_u - g_k];
    dG_k = [dg_k; -dg_k];
    % Number of constraints
    nEq = length(h(x0));
    nIneq = length(G_k);

    % Storing info values
    x = zeros(n, MaxIter);
    df= zeros(n, MaxIter);
    f = zeros(1, MaxIter);
    % Equality constraint
    y = zeros(nEq, MaxIter);
    % Inequality constraints
    z = zeros(nIneq, MaxIter);
    % Box constraints
    a = zeros(n, MaxIter);
    b = zeros(n, MaxIter);

    % Setting Initial values
    x(:,1) = x0;
    [f(1), df(:,1)] = objective(x0);
    
    % Initial Hessian
    H = zeros(n, n, MaxIter); 
    if nargin < 10
        H(:,:,1) = 10*eye(n); % Default to identity matrix 
    else
        H(:,:,1) = hessianfunc(x0, y(:, 1), z(:, 1));
    end

    fprintf('initial Hessian\n');
    disp(H(:,:,1));

    % Main loop

    while k < MaxIter && ~converged

        k = k+1;

        fprintf('-------- Iteration %d --------\n', k);
        fprintf('Current x: \n');
        disp(x(:,k));
        %{
            Solving local QP for search direction p_k:

            min 1/2 p_k' H(x_k) p_k + ∇ f(x_k)' p_k + f(x_k)

            s.t:
                ∇ h(x_k)' p_k + h_k = 0
                ∇ g(x_k)' p_k + g_k = 0

            quadprog requires:
                min 1/2 p_k' H p_k + grad^T p_k
                s.t:
                    Aeq p_k = beq
                    Aineq p_k <= bineq
                    lb <= p_k <= ub

            [x,fval,exitflag,output,lambda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options) 

        %}
        % Linearised constraints
        
        % A = -dG_k'; THIS SHOULD BE TRANSPOSED, NOT SURE WHY I AM GETTING DIMENSION ERRORS
        % PROCEEDING ANYWAY
        A  = -dG_k;
        b_quadprog = G_k; % Had to be explicitly defined as b_quadprog to avoid conflict with b
        Aeq = dh_k';
        beq = -h_k;

        % fprintf('x\n');
        % disp(x(:,k));
        % fprintf('x_l\n');
        % disp(x_l);

        % Lower and upper bounds for search direction 
        p_k_ub = x_u - x(:,k);
        p_k_lb = x_l - x(:,k);


        fprintf('-------- Solving QP %d------\n', k);
        fprintf('Hessian for QP: \n');
        disp(H(:,:,k));

        [p_k, ~, ~, ~, l] = quadprog(H(:,:,k), df(:,k), A, b_quadprog, Aeq, beq, p_k_lb, p_k_ub);
        
        fprintf('-------- QP Complete --------\n');

        fprintf('Search direction: \n');
        disp(p_k);

        % Defining Lagrange multiplier directions
        dy = l.eqlin - y(:,k);
        dz = l.ineqlin - z(:,k);
        da = l.lower - a(:,k);
        db = l.upper - b(:,k);

        % Line search
        if strcmp(method, 'line-search')
            alpha = lineSearch(x(:,k), objective, h, g, g_l, g_u, p_k, y(:,k), z(:,k));

        elseif strcmp(method, 'trust-region')
            % Trust region
            alpha = trustRegion(x(:,k), objective, h,g, g_l, g_u , y(:,k), z(:,k), dx);
        elseif strcmp(method, 'local')
            % Local method
            alpha = 1;
        else
            error('Invalid method. Please choose either line-search or trust-region');
        end

        fprintf('Alpha: %f\n', alpha);

        % Updating lagrange multipliers
        y_new = y(:,k) + alpha*dy;
        z_new = z(:,k) + alpha*dz;
        a_new = a(:,k) + alpha*da;
        b_new = b(:,k) + alpha*db;

        % z = inequality 
        % y = equality

        % ----------- Final update steps --------------- (Slide 21, Lec SQP)
        % Compute: ∇x L(xk , yk+1, zk+1) = ∇f (xk ) − ∇h(xk )yk+1 − ∇g(xk )zk+1

        
        D_x_Lagrange_1 = df(:,k) - dh_k*y_new - dG_k'*z_new - a_new + b_new;

        fprintf('D_x_Lagrange_1\n');
        disp(D_x_Lagrange_1);


        %Compute: xk+1 = xk + ∆xk
        
        x(:,k+1) = x(:,k) + alpha*p_k;
        
        % Update: f (xk+1), ∇f (xk+1), h(xk+1), ∇h(xk+1), g(xk+1 ),∇g(xk+1)
        
        [f(k+1), df(:,k+1)] = objective(x(:,k+1));
        [h_k, dh_k] = h(x(:,k+1));
        [g_k, dg_k] = g(x(:,k+1));
        G_k = [g_k - g_l; g_u - g_k];
        dG_k = [dg_k; -dg_k];

        % Compute: ∇xL(xk+1, yk+1, zk+1) = ∇f (xk+1) − ∇h(xk+1)yk+1 − ∇g(xk+1)zk+1

    
        
        D_x_Lagrange_2 = df(:,k+1) - dh_k*y_new - dG_k'*z_new - a_new + b_new;

        fprintf('D_x_Lagrange_2\n');
        disp(D_x_Lagrange_2);
        
        % Compute pk = xk+1 − xk and qk = ∇x L(xk+1, yk+1, zk+1) − ∇x L(xk , yk+1, zk+1)
        
        p = x(:,k+1) - x(:,k);
        q = D_x_Lagrange_2 - D_x_Lagrange_1;

        fprintf('------ Checking dimensions of p and q ------\n');
        fprintf('p\n');
        disp(p);
        fprintf('q\n');
        disp(q);

        fprintf('------ Performing BFGS update ------\n');
        
        % Update Hessian
        
        if nargin < 10
            H(:,:,k+1) = dampedBFGS(H(:,:,k), p, q);
            fprintf('---- Updated Hessian: ------- \n');
            disp(H(:,:,k+1));
        else
            H(:,:,k+1) = hessianfunc(x(:,k+1), y_new, z_new);
        end

        k  = k + 1;

        converged = norm(D_x_Lagrange_2, inf) < tol && ...
                    norm(p, inf) < tol && ...
                    norm(q, inf) < tol && ...
                    norm(dy, inf) < tol && ...
                    norm(dz, inf) < tol && ...
                    norm(da, inf) < tol && ...
                    norm(db, inf) < tol;

    end

    % Output
    fval = f(k);

    y(:,k) = y_new;
    z(:,k) = z_new;
    a(:,k) = a_new;
    b(:,k) = b_new;

    lambda.eq = y(:,k);
    lambda.ineq = z(:,k);
    lambda.lower = a(:,k);
    lambda.upper = b(:,k);

    grad = df(:,k);
    H = H(:,:,k);

    iterations = k;

    % Output message

    fprintf('-------- SQP Complete ---------\n');
    fprintf('Optimal solution: \n');
    disp(x(:,k));
    fprintf('Optimal function value: \n');
    disp(fval);
    fprintf('Number of iterations: \n');
    disp(iterations);
    fprintf('Lagrange multipliers: \n');
    fprintf('Equality constraints: \n');
    disp(lambda.eq);
    fprintf('Inequality constraints: \n');
    disp(lambda.ineq);
    fprintf('Lower box constraints: \n');
    disp(lambda.lower);
    fprintf('Upper box constraints: \n');
    disp(lambda.upper);
    fprintf('Final gradient: \n');
    disp(grad);
    fprintf('Final Hessian: \n');
    disp(H);

    if ~converged
        exitflag = 0;
        fprintf('Maximum number of iterations reached, without convergence.\n');
    else
        exitflag = 1;
        fprintf('Succesfully converged.\n');
    end

end
