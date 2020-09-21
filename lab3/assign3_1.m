clear all
clc
clf
close all
x = ones(2, 1);
solve(x, 1);
x = zeros(2, 1);
solve(x, 2);
function solve(x, t)
    del = ones(2, 1);
    iter=0;
    fprintf('iter\t\tx\t\ty\t\tf1\tf2\n');
    while(iter<2)
        j = jacobian(x(1), x(2), t);
        f = calculate(x(1), x(2), t);
        fprintf('%d\t%14d \t%14d\t%14d\t%14d\n', iter, x(1), x(2), f(1), f(2));
        del = mldivide(j, f);
        x = x+del;
        iter = iter+1;
    end
    f = calculate(x(1), x(2), t);
    fprintf('%d\t%14d \t%14d\t%14d\t%14d\n', iter, x(1), x(2), f(1), f(2));
    fprintf('\n\n');
end
 
function f = calculate(x, y, t)
    f=zeros(2, 1);
    if(t==1)
        f(1) = x*y*y + x*x*y + x^4-3; 
        f(2) = (x^3)*(y^5)-2*(x^5)*y - x*x+2;
    end
    if(t==2)
        f(1) = sin(4*pi*x*y)-2*y-x;
        f(2) = ((4*pi-1)/4*pi)*(exp(2*x)-exp(1))-2*exp(1)*x+4*(exp(1))*y*y;
    end
end
function j = jacobian(x, y, t)
    j = zeros(2, 2);
    if(t==1)
        j(1, 1) = y*y + 2*x*y + 4*x*x*x;
        j(1, 2) = 3*x*x*(y^5) - 10*(x^4)*y-2*x;
        j(2, 1) = 2*x*y+x*x;
        j(2, 2) = 5*(x^3)*(y^4)-2*(x^5);
    end
    if(t==2)
        j(1, 1) = cos(4*pi*x*y)*4*pi*y - 1;
        j(1, 2) = ((4*pi-1)/4*pi)*(2*exp(2*x))-2*exp(1);
        j(2, 1) = cos(4*pi*x*y)*4*pi*x - 2;
        j(2, 2) = 8*exp(1)*y;
    end
end