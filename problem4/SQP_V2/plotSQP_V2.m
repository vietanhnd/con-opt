%{
    Script for testing SQP algorithm.

    Compares the results of our SQP algorithm with the results from fmincon.

    1. Himmeblau's test problem

%}

% Test SQP on Himmelblau's test problem

% Initial guess
x0 = [0; 0];
% Lower and upper bounds
x_l = [-5; -5];
x_u = [5; 5]; 
% Lower and upper bounds for inequality constraints
g_l = 0;
g_u =  1e6; % Large number for upper bound (inf)
% Equality constraint

% h = @himmelblauEq;
h = @noEqCon;
% Inequality constraint
g = @himmelblauIneq;
% Objective function
objective = @himmelblau;
% % Hessian function - USING BFGS
% hessianfunc = @himmelblauHessian;
% % Lagrange Hessian function
HBlagrangeHessian = @himmelblauLagrangeHessian;

% Run SQP
method = 'line-search';
% x, fval, exitflag, iterations, lambda, grad, H] = SQP(objective, x0, x_l, x_u, g_l, g_u, h, g, method);
[x, fval, exitflag, iterations, lambda, grad, H] = SQP_V2(objective, x0, x_l, x_u, g_l, g_u, h, g, method);

% Plotting results


plotHBwIter(x, false);


% print results
function [f, df] = himmelblau(x)
    f = (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
    df = [4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 - 7); 2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7)];
end

% function H = himmelblauLagrangeHessian(x, lambdaEq, lambdaIneq)

% end
function [h, dh, d2h] =  himmelblauEq(x)
    %{
        Equality constraint for Himmelblau's test problem.

        y = 2/3 x
    %}

    h = [(2/3)*x(1) - x(2)];
    dh = [2/3;-1];

    % REMOVING EQUALITY CONSTRAINT FOR NOW

    % h = 0;
    % dh = [0;0];
    d2h = [[0 0;0 0], [0 0;0 0]];
end

function [h, dh, d2h] = noEqCon(x)
    h = 0;
    dh = [0;0];
    d2h = [[0 0;0 0], [0 0;0 0]];
end


function [g, dg, d2g] = himmelblauIneq(x)
    g = [(x(1)+2)^2 - x(2); -4*x(1) + 10*x(2)];
    dg = [2*(x(1)+2), -4; -1, 10];
    d2g = [[2 0;0 0], [0 0;0 0], [-2 0;0 0], [0 0;0 0]];
end

function H = himmelblauHessian(x)
    %{
        Hessian of the objective function for our Himmelblau test problem.
    %}

    H = [12*x(1)^2 + 4*x(2) - 42, 4*x(1) + 4*x(2); 4*x(1) + 4*x(2), 4*x(1) + 12*x(2)^2 - 26];
end


function d2L = himmelblauLagrangeHessian(x, y, z)
    %{
        Hessian of the Lagrangian for our Himmelblau test problem.
    %}

    G_1 = [2 0; 0 0];
    G_3 = [-2 0;0 0];
    d2L = himmelblauHessian(x) - z(1)*G_1 - z(3)*G_3;

end




function plotHBwIter(x_iter, eqcon)

    % Plotting iterations on himmelblau contour plot along with feasible region and equality constraints
    % x: data from iterations


    x1 = linspace(-6, 6, 100);
    x2 = linspace(-6, 6, 100);
    [X1, X2] = meshgrid(x1, x2);


    Z = zeros(size(X1));
    for i = 1:size(X1, 1)
        for j = 1:size(X1, 2)
            Z(i, j) = himmelblau([X1(i, j); X2(i, j)]);
        end
    end

    % Plot the contour
    figure;
    v = [-10:0.5:0, 0:1:10 10:5:100 100:20:200, 200:50:400];
    contour(X1, X2, Z, v);

    % Shading the feasible region

    yc1 = (x1+2).^2;
    yc2 = 0.4.*x1;

    % Adding equality constraint: y = x, to plot

    % Adding constraints+feasible regions to the plot

    hold on

        % Equality constraint that goes through the feasible minimum, but not through the origin

        if eqcon 
            yeq = 2/3*x1;

            plot(x1, yeq, 'LineWidth', 1, 'Color', 'black');
        end
    
        
        % Inequality constraints
        fill(x1,yc1,[0.7 0.7 0.7],'facealpha',0.5)
        fill([x1 x1(end) x1(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.5)

        % Box Constraints
        % Shading area below lower bound
        fill([x1(1), x1(1), x1(end), x1(end)], [-6, -5, -5, -6], [0.7 0.7 0.7], 'facealpha', 0.5)
        % Shading area above upper bound
        fill([x1(1), x1(1), x1(end), x1(end)], [5, 6, 6, 5], [0.7 0.7 0.7], 'facealpha', 0.5)
        % Shading area to the left of the left bound
        fill([-6, -5, -5, -6], [x2(1), x2(1), x2(end), x2(end)], [0.7 0.7 0.7], 'facealpha', 0.5)
        % Shading area to the right of the right bound
        fill([5, 6, 6, 5], [x2(1), x2(1), x2(end), x2(end)], [0.7 0.7 0.7], 'facealpha', 0.5)
        
    hold off

    xlim([-6 6])
    ylim([-6 6])
    colorbar
    xlabel('x_1');
    ylabel('x_2');

    hold on

    x = x_iter(1, :);
    y = x_iter(2, :);

    plot(x, y, '-x', 'Color', [0.8500, 0.3250, 0.0980], 'MarkerSize', 10, 'LineWidth', 2)

    % Mark final point
    % Blue RBG: [0, 0.4470, 0.7410]
    plot(x(end), y(end), '+', 'Color', [0, 0.4470, 0.7410], 'MarkerSize', 10, 'LineWidth', 2)



end