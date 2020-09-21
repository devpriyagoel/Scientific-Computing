clear all;
clc;
clf;
close all;

fprintf('\n')
X1 = solve(-1, 1.5);
X2 = solve(-1, 1.7);
X3 = solve(-1, 1.2);
function X = solve(x, c)
    i=0;
    fprintf('\nFor c = %f\n',c);
    fprintf('\tn\tx(n)\t\t\tf(x(n))\n\n')
    while(abs(f(x,c))>1e-6)
        fprintf('\t%d\t%e\t\t%13e \t\n',i,x,f(x,c));
        X(i+1,1)=i;
        X(i+1,2)=x;
        X(i+1,3)=f(x,c);
        xx = x - f(x,c)/f_dash(x);
        x = xx;
        i=i+1;
    end
    fprintf('\n')
end

function y = f(x,c)
    y = exp(x)-c-atan(x);
end

function y = f_dash(x)
    y = exp(x)-(1/(1+x^2));
end