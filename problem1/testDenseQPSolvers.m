
%{

    Initial testing of dense methods,
    
    (example 16.2, Nocedal and Wright)

    Solution:
    - x* = (2, -1, 1)'
    - λ* = (3, -2)'

%}

x_sol = [2; -1; 1];
lambda_sol = [3; -2];

% Define problem

H = [[6, 2, 1];
     [2, 5, 2];
     [1, 2, 4]];

g = [-8; -3; -3];

A = [[1, 0];
     [0, 1];
     [1, 1]];

b = [3; 0];

% Solve using dense LDL

[x, lambda] = EqualityQPSolver(H, g, A, b, 'LDLdense');

disp('LDL dense solution:');
disp(x);
disp(lambda);

% Solve using dense LU

[x, lambda] = EqualityQPSolver(H, g, A, b, 'LUdense');
disp('LU dense solution:');
disp(x);
disp(lambda);

% Solve using Null-space

[x, lambda] = EqualityQPSolver(H, g, A, b, 'nullspace');
disp('Null-space solution:');
disp(x);
disp(lambda);

% Solve using Range-space

[x, lambda] = EqualityQPSolver(H, g, A, b, 'rangespace');
disp('Range-space solution:');
disp(x);
disp(lambda);


% Testing sparse solvers

[x, lambda] = EqualityQPSolver(H, g, A, b, 'LDLsparse');
disp('LDL sparse solution:');
disp(x);
disp(lambda);

[x, lambda] = EqualityQPSolver(H, g, A, b, 'LUdense');
disp('LU sparse solution:');
disp(x);  
disp(lambda);



