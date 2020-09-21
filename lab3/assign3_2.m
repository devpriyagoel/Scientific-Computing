clear all
clc
clf
close all
%f= [@(x1, x2, x3) 6*x1 -2*cos(x2* x3 )- 1, @(x1, x2, x3) 9*x2 + sqrt(x1^2 + sin(x3) + 1.06)+0.9 , @(x1, x2, x3) 60*x3 + 3*exp( -x1* x2) + 10*pi-3];
syms x1 x2 x3
f1 = @(x1,x2,x3)6*x1-2*cos(x2*x3 )-1;
f2 = @(x1,x2,x3)9*x2+sqrt(x1^2+sin(x3)+1.06)+0.9;
f3 = @(x1,x2,x3)60*x3+3*exp(-x1*x2)+10*pi-3;

syms x1 x2 x3
Jac = jacobian([f1(x1,x2,x3),f2(x1,x2,x3),f3(x1,x2,x3)],[x1,x2,x3])
t = zeros(3, 1);
solve(t);
function solve(t)
    del = ones(3, 1);
    iter=0;
    fprintf('iter\t\tx   \t     y \t\t\tz \t\t f1 \t\tf2\t\tf3\n');
    while(norm(del, Inf)>1e-6)
        j = calc_jacobian(t(1), t(2), t(3));
        f = calc_func(t(1), t(2), t(3));
        fprintf('%d\t%14d \t%14d \t%14d \t%14d \t%14d \t%14d\n', iter, t(1), t(2), t(3), f(1), f(2), f(3));
        del = mldivide(j, f);
        t = t-del;
        iter = iter+1;
        
    end
    f = calc_func(t(1), t(2), t(3));
    fprintf('%d\t%14d \t%14d \t%14d \t%14d \t%14d \t%14d\n', iter, t(1), t(2), t(3), f(1), f(2), f(3));
    fprintf('\n Solution : %14d \t%14d \t%14d\n', t(1), t(2), t(3));
end

function y = calc_func(x1, x2, x3)
    y = zeros(3, 1);
    y(1) = 6*x1-2*cos(x2*x3 )-1;
    y(2) = 9*x2+sqrt(x1^2+sin(x3)+1.06)+0.9;
    y(3) = 60*x3+3*exp(-x1*x2)+10*pi-3;
end

function y = calc_jacobian(x1, x2, x3)
    y(1, 1) = 6;
    y(1, 2) = 2*x3*sin(x2*x3);
    y(1, 3) = 2*x2*sin(x2*x3);
    y(2, 1) = x1/(x1^2 + sin(x3) + 53/50)^(1/2);
    y(2, 2) = 9;
    y(2, 3) = cos(x3)/(2*(x1^2 + sin(x3) + 53/50)^(1/2));
    y(3, 1) = -3*x2*exp(-x1*x2);
    y(3, 2) = -3*x1*exp(-x1*x2);
    y(3, 3) = 60;
end
