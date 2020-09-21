clear all
clc
clf
close all

%x = 1+1/(x+1)x
figure(1)
fplot(@(x) 1+1./((x+1).*x), [0.5, 2]);
hold on;
fplot(@(x) x, [0.5, 2]);
x=1; y=0;
fprintf('iteration \t x \t f(x)\n');
iter=0;
while(int32(x*100) ~= int32(y*100))
    y = x;
    x = calculate(x);
    iter = iter + 1;
    line([y, y], [y, x], 'LineStyle', '-', 'Color', 'k'); 
    line([y, x], [x, x], 'LineStyle', '-', 'Color', 'k'); 
    fprintf('%5d \t %12d \t %12d\n', iter, y, x);
end
fprintf('\nSolution  = %d\n', x);
function f = calculate(x)
    f = 1+1./((x+1).*x);
end

