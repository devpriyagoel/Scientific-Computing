clear all;
clc;
clf;
close all;

X1 = solve(1);
X2 = solve(2);
X3 = solve(3);
X4 = solve(4);


function X = solve(t) 
    x0=1;
    fprintf('\n\n\t\tn\t\t\tx(n)\t\n\n')
    for i = 1:4
        X(i,1) = i-1;
        X(i,2) = calculate(x0, t);
        fprintf('\t%13e\t\t%13e\t\n',x0,X(i, 2));
        x0 = X(i, 2);
    end
end
function f = calculate(x, t)
    f=double(0);
    if t==1
        f = x*((1+((7-(x^5))/(x^2)))^3);
    end
    if t==2
        f = x-(x^5-7)/(x^2);
    end
    if t==3
        f = x-(x^5-7)/(5*(x^4));
    end
    if t==4
        f = x-(x^5-7)/(12);
    end
    
end