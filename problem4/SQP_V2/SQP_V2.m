function [x, f, exitflag, iterations, lambda, grad, H] = SQP_V2(objective, x0, x_l, x_u, g_l, g_u, h, g, method, MaxIter, tol, hessianfunc)
    % initialisation

    if nargin < 12
        tol = 1e-3;
    end
    if nargin < 11
        MaxIter = 50;
    end
    if nargin < 10 % User has not specified hessian function
        hessianfunc = @dampedBFGS;
    end

    n = length(x0);
    nEq = length(h(x0));
    nIneq = 2*length(g(x0)); % 2x to account for both lower and upper bounds

    k = 0;
    converged = false;

    x = zeros(n, MaxIter);
    f = zeros(1, MaxIter);
    df = zeros(n, MaxIter);

    lambdaEq = zeros(nEq, MaxIter);
    lambdaIneq = zeros(nIneq, MaxIter);
    lambdaLower = zeros(n, MaxIter);
    lambdaUpper = zeros(n, MaxIter);

    % Check if x0 is within bounds - if not place in the middle of the box
    if any(x0 < x_l) || any(x0 > x_u)
        fprintf('Initial guess is not feasible. Placing it in the middle of the box.\n');
        x0 = (x_l + x_u)/2;
    end

    % Initial values
    x(:,1) = x0;
    [f(1), df(:,1)] = objective(x0);

    [h_k, dh_k, d2h_k] = h(x0);
    [g_k, dg_k, d2G_k] = g(x0); % Writing inequality constraints in standard form
    G_k = [g_k - g_l; g_u - g_k];
    dG_k = [dg_k, -dg_k];



    % Initial Hessian
    H = zeros(n, n, MaxIter);
    if nargin < 10
        H(:,:,1) = 5*eye(n); % Default to identity matrix 
    else
        H(:,:,1) = hessianfunc(x0, lambdaEq(:, 1), lambdaIneq(:, 1));
    end

    % Initial Trust region radius
    Delta_k = 1;
    Delta_max = 1e6; % Not sure what this should be set to!?

    

    % Main loop
    fprintf('Starting SQP algorithm\n');


    while k<MaxIter && ~converged 

        k = k + 1;
        accept_step = false; % for trust-region method (line search always accepts final step)
   
        fprintf('Iteration %d\n', k);

        %{
            Solving QP subproblem
        %}

        % fprintf('----- Solving QP subproblem -----\n');

        A = -dG_k';
        b = G_k;
        Aeq = dh_k';
        beq = -h_k;


        % fprintf('A\n')
        % disp(A);
        % fprintf('b\n')
        % disp(b);



        % Bounds for QP subproblem - Different according to method

        p_k_lb = [];
        p_k_ub = [];

        if strcmp(method, 'trust-region')
            % Trust region
            p_k_lb = -Delta_k*ones(n, 1);
            p_k_ub = Delta_k*ones(n, 1);
            
        elseif strcmp(method, 'line-search')
            % Line search
            p_k_lb = x_l - x(:, k);
            p_k_ub = x_u - x(:, k);
        elseif strcmp(method, 'local')
            % Local method - no bounds
            
            p_k_ub = inf*ones(n, 1);
            p_k_lb = -inf*ones(n, 1);
            
        else
            error('Invalid method. Please choose either line-search or trust-region');
        end


        % Hessian for QP subproblem - Hessian of the Lagrangian

        H_lagrange = lagrangeHessian(H(:,:,k), lambdaIneq(:, k), lambdaEq(:, k),d2G_k,d2h_k);
        [p_k, ~, ~, ~, l] = quadprog(H_lagrange, df(:,k), A, b, Aeq, beq, p_k_lb, p_k_ub);


        % Direction of Lagrange multiplier step
        p_lambdaEq = l.eqlin - lambdaEq(:, k);
        p_lambdaIneq = l.ineqlin - lambdaIneq(:, k);
        p_lambdaLower = l.lower - lambdaLower(:, k);
        p_lambdaUpper = l.upper - lambdaUpper(:, k);

        
        % Penalty parameter for line search and trust-region
        %{
            μ ≥ ∇fk'pk+(σ/2)pk'∇2xxLk pk/(1 − ρ)‖ck ‖1

        %}
        % Checking if lagrange hessian is positive definite
        sigma = 1;
        if any(eig(H_lagrange) <= 0)
            % Lagrange Hessian not positive definite
            sigma = 0;
        end

        mu = (df(:, k)'*p_k + (sigma/2)*p_k'*H_lagrange*p_k)/(1 - 0.5*norm([h_k;G_k], 1));



        fprintf('Penalty parameter: %f\n', mu);
        alpha = 1;
        

        if strcmp(method, 'line-search')
            %{
                Back-tracking Line search

                pλ = ˆλ − λk Choose μk ≥ ∇fk pk +(1/2)p′k Bk pk /(1−ρ)‖ck‖1 and set αk ← 1 
                
                while φ1(xk + αk pk ; μk ) > φi(xk ; μk ) + ηαk D(φ1((xk ; μk ); pk )) do 
                    αk ← τα αk for some τα ∈ (0, τ ] 
                end while
            %}

            fprintf('----- Performing Line search -----\n');

   
            eta = 0.25;
            tau = 0.5;

            % φ1 (x; μ) = f (x) + μ‖c(x)‖1. 
            % D(φ1 (x k ; μ); p k ) = ∇ f T k p k − μ‖ck ‖1.
            
            x_test = x(:, k) + alpha*p_k;
            [f_test, ~] = objective(x_test);
            [h_test, ~, ~] = h(x_test);
            [g_test, ~, ~] = g(x_test);
            G_test = [g_test - g_l; g_u - g_test];
    

            while (f_test + mu*norm([h_test;G_test], 1)) > (f(k)+ mu*norm([h_test;G_test], 1))+eta*alpha*(df(:,k)'*p_k-mu*norm([h_test;G_test], 1))
                alpha = tau*alpha;
                x_test = x(:, k) + alpha*p_k;
                [f_test, ~] = objective(x_test);
                [h_test, ~, ~] = h(x_test);
                [g_test, ~, ~] = g(x_test);
                G_test = [g_test - g_l; g_u - g_test];
            end
            accept_step = true;

            fprintf('---- Line Search Complete ----\n')
            fprintf('Alpha: %f\n', alpha);
    
        elseif strcmp(method, 'trust-region')
            %{
                Trust region
            %}
            fprintf('----- Performing Trust region -----\n');
            eta = 0.1; % Not sure about this?!
         
            % q_u_0 = f(k) + mu*norm(h_k, 1) + mu*sum(max(0, -(G_k + dG_k'*p_k)));
            

         
            
            % % fprintf('q_u_0: %f\n', q_u_0);

        
            % phi_x = f(k) + mu*norm(h_k, 1)+mu*sum(max(0, -(G_k + dG_k'*p_k)));
            
            % % Getting updated function value and constraints, according to step pk

            % x_pk = x(:, k) + p_k;
            % [f_pk, ~] = objective(x_pk);
            % [h_pk, ~, ~] = h(x_pk);
            % [g_pk, ~, ~] = g(x_pk);
            % G_pk = [g_pk - g_l; g_u - g_pk];

            % phi_x_pk = f_pk + mu*norm(h_pk, 1)+mu*sum(max(0, -G_pk));

            % q_u_pk = f_pk + mu*norm(h_pk, 1) + mu*sum(max(0, -(G_pk + dG_k'*p_k)));

            % rho = (phi_x - phi_x_pk)/(q_u_0 - q_u_pk);

            % % Determining how to update trust region radius
            % % fprintf('Rho: %f\n', rho);
            % % fprintf('||p_k||_inf: %f\n', norm(p_k, inf));

            % if rho < 0.25
            %     Delta_k = 0.25*Delta_k;
            % else
            %     if rho > 0.75 && norm(p_k, inf) == Delta_k
            %         Delta_k = min(2*Delta_k, Delta_max);
            %     else
            %         Delta_k = Delta_k;
            %     end
            % end
            
            % % Determining if step is accepted
            % if rho > eta
            %     fprintf('Step accepted\n');
            %     accept_step = true;
            % end


        elseif strcmp(method, 'local')
            %{
                Local method
            %}
            fprintf('----- Performing Local method (Accepting Step) -----\n');
            accept_step = true;

        else 
            error('Invalid method. Please choose either line-search or trust-region');
        end

        if accept_step

            x(:, k+1) = x(:, k) + alpha*p_k;
            % Update values
            fprintf('---- Updating values ----\n');
            fprintf('New x: \n');
            disp(x(:, k+1));

            % Updating lagrange multipliers

            lambdaEq(:, k+1) = lambdaEq(:, k) + alpha*p_lambdaEq;
            lambdaIneq(:, k+1) = lambdaIneq(:, k) + alpha*p_lambdaIneq;
            lambdaLower(:, k+1) = lambdaLower(:, k) + alpha*p_lambdaLower;
            lambdaUpper(:, k+1) = lambdaUpper(:, k) + alpha*p_lambdaUpper;

            % ∇xL(xk , λk+1) = ∇f (xk ) − ∑ i ∇ci(xk )(λk+1)i
            % Update of Lagrange w.r.t. Lagrange multipliers
    

            D_x_Lagrange_1 = df(:, k) -  dh_k*lambdaEq(:, k+1) - dG_k*lambdaIneq(:, k+1) - lambdaLower(:, k+1) + lambdaUpper(:, k+1);

            [f(k+1), df(:, k+1)] = objective(x(:, k+1));
            [h_k, dh_k, d2h_k] = h(x(:, k+1));
            [g_k, dg_k, d2G_k] = g(x(:, k+1)); % Writing inequality constraints in standard form
            G_k = [g_k - g_l; g_u - g_k];
            dG_k = [dg_k, -dg_k];

            % ∇xL(xk+1, λk+1) = ∇f (xk+1) − ∑ i ∇ci(xk+1)(λk+1)i
            % Update of Lagrange w.r.t. Lagrange multipliers and x

            D_x_Lagrange_2 = df(:, k+1) -  dh_k*lambdaEq(:, k+1) - dG_k*lambdaIneq(:, k+1) - lambdaLower(:, k+1) + lambdaUpper(:, k+1);

            % Compute pk = xk+1 − xk and qk = ∇x L(xk+1, yk+1, zk+1) − ∇x L(xk , yk+1, zk+1)
            p = x(:, k+1) - x(:, k); % Should be equal to alpha*p_k
            q = D_x_Lagrange_2 - D_x_Lagrange_1;

            % Update Hessian
            if nargin < 10
                % BFGS update
                B = H(:,:,k);
                theta = 1;
                if p'*q < 0.2*p'*B*p
                    theta = 0.8*p'*(B*p)/(p'*(B*p) - p'*q);
                end
                
                r = theta*q + (1-theta)*B*p;

                Bk = B + (B*p)*(p'*B)/(p'*B*p) + (r*r')/(p'*r);
                H(:,:,k+1) = Bk;
            else
                H(:,:,k+1) = hessianfunc(x(:, k+1), lambdaEq(:, k+1), lambdaIneq(:, k+1));
            end


            fprintf('---- Updated Hessian ---\n');
            fprintf('Previous Hessian\n');
            disp(H(:,:,k));
            fprintf('Updated Hessian\n');
            disp(H(:,:,k+1));

            converged = norm(D_x_Lagrange_2, inf) < tol && ...
                norm(p, inf) < tol && ...
                norm(q, inf) < tol && ...
                norm(dy, inf) < tol && ...
                norm(dz, inf) < tol && ...
                norm(da, inf) < tol && ...
                norm(db, inf) < tol;

            

        % else
        %     % Step not accepted
        %     x(:, k+1) = x(:, k);
        %     f(k+1) = f(k);
        %     df(:, k+1) = df(:, k);
        %     lambdaEq(:, k+1) = lambdaEq(:, k);
        %     lambdaIneq(:, k+1) = lambdaIneq(:, k);
        %     lambdaLower(:, k+1) = lambdaLower(:, k);
        %     lambdaUpper(:, k+1) = lambdaUpper(:, k);
        %     H(:,:,k+1) = H(:,:,k);
        end
    
        % Check convergence - Not sure if this is correct
  
        % % Dummy converged condition - remove this!
        % converged = true;
    end


    if converged
        exitflag = 1;
    else
        exitflag = 0;
    end

    iterations = k;
    lambda.eq = lambdaEq(:, k);
    lambda.ineq = lambdaIneq(:, k);
    lambda.lower = lambdaLower(:, k);
    lambda.upper = lambdaUpper(:, k);
    grad = df;

    fprintf('-------- SQP Complete ---------\n');
    if exitflag == 1
        fprintf('Program successfully converged\n');
    else
        fprintf('Program did not converge\n');
    end
    fprintf('Optimal solution: \n');
    disp(x(:,k));
    fprintf('Optimal function value: \n');
    disp(f(k));
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
    disp(grad(:,k));
    fprintf('Final Hessian: \n');
    disp(H(:,:,k));





end



function Lxx = lagrangeHessian(H, lambdaIneq, lambdaEq, d2G, d2h)
    %{
        Calculates the Hessian of the Lagrangian function.
        
        Used for QP subproblem.
    %}

    sumh = 0;
    sumG = 0;
    for i = 1:length(lambdaEq)
        sumh = sumh + lambdaEq(i)*d2h(:, i);
    end
    for i = 1:length(lambdaIneq)
        sumG = sumG + lambdaIneq(i)*d2G(:, i);
    end

    Lxx = H - sumG - sumh;

end

