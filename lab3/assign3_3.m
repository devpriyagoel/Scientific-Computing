clear all
clc
clf
close all
syms x
f = @(x)x^3-(2*x*(4-tan(x)))/(3*x*x+sin(x))+((4-tan(x))/(3*x*x+sin(x)))^7- (4*(x^3)*(4-tan(x)))/(3*x*x+sin(x))-5;
syms x
Jac = jacobian([f],[x])
calc_y = @(x) (4-tan(x))/(3*x*x + sin(x));
f1 = @(x,y) x^3 - 2*x*y + y^7 - 4*x^3*y-5;
f2 = @(x,y) y*sin(x)+3*x*x*y + tan(x)-4;
t = 1;
solve(t, f, calc_y, f1, f2);
function solve(t, func, calc_y, f1, f2)
    del = 1;
    iter=0;
    fprintf('iter\t\tx   \t     y \t\t f1(x, y) \t f2(x, y)\n');
    while(abs(del)>1e-6)
        j = diff(t);
        f = func(t);
        y = calc_y(t);
        fprintf('%d\t%13d\t13%d\t%13d\t%13d\n', iter, t, y, f1(t, y), f2(t, y));
        del = f/j;
        t = t-del;
        iter = iter+1;        
    end
    y = calc_y(t);
    fprintf('%d\t%13d\t%13d\t%13d\t%13d\n', iter, t, y, f1(t, y), f2(t, y));
    fprintf('\n Solution : %14d %14d \n', t, y);
end


function y=diff(x)
    y=(2*(tan(x) - 4))/(sin(x) + 3*x^2) + 3*x^2 + (2*x*(tan(x)^2 + 1))/(sin(x) + 3*x^2) + (12*x^2*(tan(x) - 4))/(sin(x) + 3*x^2) - (7*(tan(x) - 4)^6*(tan(x)^2 + 1))/(sin(x) + 3*x^2)^7 + (4*x^3*(tan(x)^2 + 1))/(sin(x) + 3*x^2) + (7*(6*x + cos(x))*(tan(x) - 4)^7)/(sin(x) + 3*x^2)^8 - (2*x*(6*x + cos(x))*(tan(x) - 4))/(sin(x) + 3*x^2)^2 - (4*x^3*(6*x + cos(x))*(tan(x) - 4))/(sin(x) + 3*x^2)^2;
end
