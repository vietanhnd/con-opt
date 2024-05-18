%{

Solves Himmelblaus test problem, using Matlab's fmincon function.

min f(x)=(x1^2 +x2 -11)^2 +(x1 +x2^2 -7)^2
s.t:


(x1+2)^2 -x2 >= 0, nonlinear inequality constraint

-4x1 +10x2 >= 0, linear inequality constraint

and 

-5 <= x1 <= 5
-5 <= x2 <= 5

%}


x0 = [0,0]; % Initial Guess - Midpoint of domain

algorithm = 'sqp';
fp = '/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem4/data/fmincon-iter';


disp(['Solving with ', algorithm])
[xsol,fval,exitflag,output,lambda,grad,hessian, history] = runfmincon_hb(x0, algorithm);
iterates=history.x;
% saving the results to a csv

filename = strcat(fp, '/results_', algorithm, '.csv');
writematrix(iterates, filename);

    

% Functions

function [c, ceq, gc, gceq] = himmelnonlcon(X)
    x = X(1);
    y = X(2);
    c = -((x+2)^2 - y);
    ceq = [];
    if nargout > 2
        gc = [-2*(x+2); 1];
        gceq = [];
    end
end


function Hout = hessianfcn(X, lambda)
    x = X(1);
    y = X(2);
    himmelblauHessian = [12*x^2 + 4*y - 42, 4*x + 4*y; 4*x + 4*y, 12*y^2 + 4*x - 26];
    Hg = [-2, 0; 0, 0];
    Hout = himmelblauHessian + lambda.ineqnonlin(1)*Hg;
end


function [xsol,fval,exitflag,output,lambda,grad,hessian, history] = runfmincon_hb(x0, algorithm)
    %{
    Solves Himmelblaus test problem, using Matlab's fmincon function.

    Allows the user to specify the algorithm used and initial guess, plotting and returning iterates.

    %}
    % Set up shared variables with outfun
    history.x = [];
    history.fval = [];
     
    % Defining Himmelblau's problem with constraints

    % Linear inequality constraint, Ax <= b (notice that the inequality is reversed compared to the usual notation)

    A = [4,-10]; 
    b = 0; 

    % No linear equality constraints

    Aeq = []; 
    beq = []; 

    % Lower and upper bounds
    lb = [-5, -5];
    ub = [5, 5];

    options = optimoptions('fmincon',...
    'OutputFcn',@outfun,...
    'Display', 'iter', 'Algorithm',algorithm,...
    'SpecifyObjectiveGradient', true,...
    'SpecifyConstraintGradient', true,...
    'HessianFcn', @hessianfcn); % SQP and active-set do not use Hessian (ignored)

    [xsol,fval,exitflag,output,lambda,grad,hessian] = fmincon(@himmelblau,x0,A,b,Aeq,beq,lb,ub,@himmelnonlcon, options);
     
     function stop = outfun(x,optimValues,state)
         stop = false;
     
         switch state
             case 'init'
                 hold on
             case 'iter'
             % Concatenate current point and objective function
             % value with history. x must be a row vector.
               history.fval = [history.fval; optimValues.fval];
               history.x = [history.x; x];

               plot(x(1),x(2),'o');

             % Label points with iteration number and add title.
             % Add .15 to x(1) to separate label from plotted 'o'.
               text(x(1)+.15,x(2),... 
                    num2str(optimValues.iteration));
               title(algorithm);
             case 'done'
                 hold off
             otherwise
         end
     end
end