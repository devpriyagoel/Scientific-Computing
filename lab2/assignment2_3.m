clear all
clc
clf
close all

solve(1, 1, 1);
solve(1, 2, 2);
solve(1, 3, 3);
solve(2, 1, 4);
solve(2, 2, 5);
solve(2, 3, 6);
 
function solve(x, t, p)
    y=0;
    figure(p);
    fplot(getFunction(t));
    hold on;
    fplot(getFunction(0));
    iter=0;
    fprintf('\niteration \t x \t f(x)\n');
    while(int32(x*1e5) ~= int32(y*1e5)&& iter<m=10)
        y = x;
        x = calculate(x, t);
        iter = iter+1;
        line([y, y], [y, x], 'LineStyle', '-', 'Color', 'k'); 
        line([y, x], [x, x], 'LineStyle', '-', 'Color', 'k'); 
        fprintf('%5d \t %12d \t %12d\n', iter, y, x);
    end
    fprintf('\nSolution = %d\n\n', x);
end
function f = calculate(x, t)
    f=double(0);
    if(t==1)
        f = x*x-1;
    end
    if(t==2)
        f = 1+2*x-x*x;
    end
    if(t==3)
        f = (1/2)*(1+3*x-x*x);
    end
end
function f = getFunction(t)
    f = @(x) x;
    if(t==1)
        f = @(x) x.*x-1;
    end
    if(t==2)
        f = @(x) 1+2*x-x.*x;
    end
    if(t==3)
        f = @(x) (1/2)*(1+3.*x-x.*x);
    end
end
