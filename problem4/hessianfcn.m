function Hout = hessianfcn(x, lambda)
    himmelblauHessian = [12*x(1)^2 + 4*x(2) - 42, 4*(x(1) + x(2)); 4*(x(1) + x(2)), 4*x(1) + 12*x(2)^2 - 26];
    Hg1 = [-2, 0; 0, 0];
    Hg2 = [0, 0; 0, 0];
    Hout = himmelblauHessian + lambda.ineqnonlin(1)*Hg1 + lambda.ineqnonlin(2)*Hg2;
end