clear all
clc
clf
close all

X1 = bisection_method(-1.5*pi, -pi, 0, 0.01);
X2 = bisection_method(0.5, 1.5, 1, 0.001);
X3 = bisection_method(0, 1, 2, 1e-3);
X4 = bisection_method(0.000, 1.000, 3, 1e-3);
X5_1 = bisection_method(-3.000, -2.000, 4, 1e-3);
X6_1 = bisection_method(-1.000, 0.000, 4, 1e-3);
X7_1 = bisection_method(0.200, 0.300, 5, 1e-3);
X7_2 = bisection_method(1.200, 1.300, 5, 1e-3);
X8 = bisection_method(2.9000, 2.9500, 6, 1e-4);

function X = bisection_method(x, y, t, error)
    i=0;
    fprintf('\t\tn\t\tx(n)\t\t\tf(x(n))\n')
    while (x<y-error)
        m = (x+y)/2;
        X(i+1,1)=i;
        X(i+1,2)=m;
        X(i+1,3)=calculate(m, t);
        fprintf('\t%13d\t%8.7e\t\t%13e \t\n',i,m,X(i+1, 3));
        if (X(i+1, 3)-error<=0 && X(i+1, 3) +error>0)
            break;
        elseif (X(i+1, 3)*calculate(y, t)<0)
            x = m;
        else
            y = m;
        end
        i=i+1;
    end
    fprintf('solution = %d \t iterations = %d', x, i); 
    fprintf('\n');
end


function f = calculate(x, t)
    f=double(0);
    if t==0
        f = exp(x)-sin(x);
    end
    if t==1
        f = 2 + cos(exp(x)-2) - exp(x);
    end
    if t==2
        f = x - pow2(-x);
    end
    if t==3
        f = exp(x)-x*x + 3*x -2;
    end
    if t==4
        f = 2*x*cos(2*x) - (x + 1)*(x+1);
    end
    if t==5
        f = x*cos(x)-2*x + 3*x-1;
    end
    if t==6
        f= x*x*x-25;
    end
    
end
