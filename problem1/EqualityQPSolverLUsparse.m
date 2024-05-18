%{
    EqualityQPSolverLUsparse:
    
    Solves Quadratic Program with sparse KKT matrix using slightly modified LU factorisation

        min 1/2x'Hx + gx subject to A'x = b, x ≥ 0

    Require:
    - H, g, A, b 
    Ensure:
    - Optimal solution x^* and λ^* found
%}

function [x, lambda] = EqualityQPSolverLUsparse(H,g,A,b)
    
[n, m] = size(A);
KKT = [[H, -A]; [-A', zeros(m)]];
v = - [g; b];

% Make KKT and RHS sparse
KKT = sparse(KKT);
v = sparse(v);


[L, U, P] = lu(KKT, "vector");
sol = U \ (L \ (v(P)));

lambda = sol(n+1:end);
x = sol(1:n); 
    
end

% function [x, lambda] = EqualityQPSolverLUsparse(H,g,A,b)

%     [n, m] = size(A);
%     KKT = [[H, -A]; [-A', zeros(m)]];
%     v = - [g; b];

%     % Convert KKT to sparse matrix
%     KKT = sparse(KKT);

%     % LU factorisation

%     [L,U,P,Q] = lu(KKT); % P and Q are permutation matrices - More efficient
    
%     v_permuted = P*v;

%     y = L \ v_permuted;
%     x = U \ y;

%     % Post multiplying x so that values are in the correct order

%     x = Q*x;


%     lambda = x(n+1:end);
%     x = x(1:n); 
% end