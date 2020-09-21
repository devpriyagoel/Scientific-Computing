clear all;
clc;
clf;
close all;

Xa = solve(1.5, 0);
Xb = solve(1.65, 1);
Xc_1 = solve(2.5, 2);
Xc_2 = solve(3.5, 2);
Xd_1 = solve(1.5, 3);
Xd_2 = solve(3, 3);
Xe_1 = solve(0.5, 4);
Xe_2 = solve(4, 4);
Xf_1 = solve(0.5, 5);
Xf_2 = solve(3.5, 5);
Xf_3 = solve(6.5, 5);

function X = solve(x, t)
    i=0;
    fprintf('\tn\tx(n)\t\t\tf(x(n))\n\n')
    while(abs(calculate(x, t))>1e-6)

        X(i+1,1)=i;
        X(i+1,2)=x;
        X(i+1,3)=calculate(x, t);
        fprintf('\t%d\t%e\t\t%13e \t\n',i,x,X(i+1, 3));
        xx = x - X(i+1, 3)/differentiate(x, t);
        x = xx;
        i=i+1;
    end
end

function f = calculate(x, t)
    f=double(0);
    if t==0
        f = exp(x)+2^(-x)+2*cos(x)-6;
    end
    if t==1
        f = log(x-1)+cos(x-1);
    end
    if t==2
        f = 2*x*cos(2*x)-(x-2)^2;
    end
    if t==3
        f = (x-2)^2-log(x);
    end
    if t==4
        f = exp(x)-3*x*x;
    end
    if t==5
        f = sin(x)-exp(-x);
    end

    
end
function f = differentiate(x, t)
    f=double(0);
    if t==0
        f = exp(x)-log(2)*2^(-x)-2*sin(x);
    end
    if t==1
        f = 1/(x-1)-sin(x-1);
    end
    if t==2
        f = -2*(2*x*sin(2*x)-cos(2*x)+x-2);
    end
    if t==3
        f = 2*(x-2)-1/x;
    end
    if t==4
        f = exp(x)-6*x;
    end
    if t==5
        f = cos(x)+exp(-x);
    end

end