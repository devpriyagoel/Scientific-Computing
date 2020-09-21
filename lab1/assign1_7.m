clear all;
clc;
clf;
close all;
N1 = newton(-0.5);
N2 = newton(0.7);
S1 = secant(0, -1);
S2 = secant(2, 1);
function X = newton(x)
    i=0;
    fprintf('\tn\tx(n)\t\t\tf(x(n))\n\n')
    while(abs(a(x))>1e-7)
        fprintf('\t%d\t%e\t\t%13e \t\n',i,x,a(x));
        X(i+1,1)=i;
        X(i+1,2)=x;
        X(i+1,3)=a(x);
        xx = x - a(x)/a_dash(x);
        x = xx;
        i=i+1;
    end
    fprintf('\n');
end
function X = secant(x0, x1)
    i=0;
    fprintf('\tn\tx(n)\t\t\tf(x(n))\n\n')
    fprintf('\t%d\t%e\t\t%13e \t\n',i,x0,a(x0));
    X(i+1,1)=i;
    X(i+1,2)=x0;
    X(i+1,3)=a(x0);
    i=i+1;
    while(abs(a(x1))>1e-7)
        fprintf('\t%d\t%e\t\t%13e \t\n',i,x1,a(x1));
        X(i+1,1)=i;
        X(i+1,2)=x1;
        X(i+1,3)=a(x1);
        xx = x1 - a(x1)*(x1-x0)/(a(x1)-a(x0));
        x0 = x1;
        x1 = xx;
        i=i+1;
    end
    X(i+1,1)=i;
    X(i+1,2)=x1;
    X(i+1,3)=a(x1);
    i=i+1;
    fprintf('\t%d\t%e\t\t%13e \t\n',i,x1,a(x1)); 
end

function f = a(x)
    f = 230*x^4 + 18*x^3 + 9*x^2 - 221*x - 9;
end

function f = a_dash(x)
    f = 920*x^4 + 54*x^2 + 18*x - 221;
end
