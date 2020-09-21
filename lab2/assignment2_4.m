clear all
clc
clf
close all

solve(-0.8, 1);
solve(0.8, 2);


function solve(x, t)
    figure(t);
    fplot(@(x) x*x+x-0.75, [-2, 1]);
    hold on;
    fplot(@(x) x, [-2, 1]);
    fprintf('\niteration \t x \t f(x)\n');
    y=0;
    iter=0;
    while(int32(x*1e5) ~= int32(y*1e5)&& iter<50)
        y = x;
        x = calculate(x);
        iter = iter+1;
        line([y, y], [y, x], 'LineStyle', '-', 'Color', 'k'); 
        line([y, x], [x, x], 'LineStyle', '-', 'Color', 'k'); 
        fprintf('%5d \t %12d \t %12d\n', iter, y, x);
    end
    fprintf('\nSolution = %d\n\n', x);
end
function f = calculate(x)
    f=x*x+x-0.75;
end
