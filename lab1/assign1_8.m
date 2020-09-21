clear all;
clc;
clf;
close all;



X1 = secant(0,3*pi/2, 1);
X2 = secant(-1.5,1.5, 2);
X3_1 = secant(-1,0.1145, 3);
X3_2 = secant(0.1145,1, 3);
X4 = secant(-1,1, 4);
X5 = secant(-1,1, 5);

function Z = secant(x0,x1, t)
    fprintf('\n');
    fprintf('t=1 %d \n', t);
    i=0;
    fprintf('\tn\tx(n)\t\t\tf(x(n))\n\n')

    Z(i+1,1)=i;
    Z(i+1,2)=x0;
    Z(i+1,3)=calculate(x0, t);
    fprintf('\t%d\t%e\t\t%13e \t\n',i,x0,Z(i+1, 3));
    i=i+1;
    while(abs(calculate(x1, t))>1e-4)
        Z(i+1,1)=i;
        Z(i+1,2)=x1;
        Z(i+1,3)=calculate(x1, t);
        fprintf('\t%d\t%e\t\t%13e \t\n',i,x1,Z(i+1, 3));
        xx = x1 - Z(i+1, 3)*(x1-x0)/(Z(i+1, 3)-calculate(x0, t));
        x0 = x1;
        x1 = xx;
        i=i+1;
    end
    Z(i+1,1)=i;
    Z(i+1,2)=x1;
    Z(i+1,3)=calculate(x1, t);
    fprintf('\t%d\t%e\t\t%13e \t\n',i,x1,Z(i+1, 3));
    i=i+1;
end

function f = calculate(x, t)
    f=double(0);
    if t==1
        f =sin(x/2)-1;
    end
    if t==2
        f = exp(x)-tan(x);
    end
    if t==3
        f = x^3-12*x^2+3*x+1;
    end
    if t==4
        f = x^3 + 4.001*x^2 + 4.002*x + 1.101;
    end
    if t==5
        f = x^6 -x^4 +2*x^3 -3*x^2 +x-4;
    end
end
