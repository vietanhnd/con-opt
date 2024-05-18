%{
    Damped BFGS update for Hessian approximation

    Requires:
    - B: Symmetric and positive definite previous hessian approximation
    - p: change in x
    - q: change in lagrange gradient 

%}
function Bk = dampedBFGS(B, p, q)

    fprintf('-------- Damped BFGS Update --------\n');
    fprintf('Previous Hessian: %d x %d\n', size(B, 1), size(B, 2));
    disp(B);

    fprintf('p: %d x 1\n', size(p, 1));
    disp(p);

    fprintf('q: %d x 1\n', size(q, 1));
    disp(q);

    theta = 1;
    if p'*q < 0.2*p'*B*p
        theta = 0.8*p'*(B*p)/(p'*(B*p) - p'*q);
    end

    r = theta*q + (1-theta)*B*p;

    Bk = B + (B*p)*(p'*B)/(p'*B*p) + (r*r')/(p'*r);

    fprintf('Theta: %d\n', theta);


    fprintf('-------- New Hessian --------\n');
    fprintf('Bk: %d x %d\n', size(Bk, 1), size(Bk, 2));
    disp(Bk);

    



    % fprintf('-------- Regular BFGS Update (for now) --------\n');

    % Bk = B + (q*q')/(p'*q) - ((B*p)*(B*p)')/((p')*(B*p));


end