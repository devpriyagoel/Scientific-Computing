clear all
clc
clf
close all
figure(3);
fplot(@(x) (3*exp(x)*(x-1))/(9*x*x), [0.6, 5]);
figure(4);
fplot(@(x) -sin(x));
solve(1, 0, 1);
solve(0, 1, 2);
function solve(x, y, t)
    figure(t);
    fplot(getFunction(t));
    hold on;
    fplot(getFunction(0));
    iter=0;
    fprintf('\niteration \t x \t f(x)\n');
    while(int32(x*1e5) ~= int32(y*1e5))
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
    f=0;
    if(t==1)
        f = exp(x)/(3*x);
    end
    if(t==2)
        f = cos(x);
    end
end
function f = getFunction(t)
    f = @(x) x;
    if(t==1)
        f = @(x) exp(x)./(3*x);
    end
    if(t==2)
        f = @(x) cos(x);
    end
end
