clear all
clc
clf
close all
muller(1.4, 1.2, 1, 1)
muller(1.3, 1.4, 1.5, 2);
muller(0.9, 1, 1.1, 3);
muller(-2.2, -2, -1.9, 3);
muller(0.7, 0.8, 0.9, 4);
muller(1i, 1i+0.1, 0.1, 5);
muller(1i, 1i+0.1, 0.1, 6);
%fplot(@(x) x-exp(x), [-5, 2]);
function muller(xkm1, xkm2, xkm3, t)
   % y = calculate(xkm1, t) + (x-xkm1)*func2(xkm1, xkm2, t) 
    %+ (x-xkm1)*(x-xkm2)*func3(xkm1, xkm2, xkm3, t); 
    iter=0;
    fprintf('\niteration \t\t x \t\t\t  f(x)\n');
    while(abs(xkm2-xkm1)>1e-5)
        w = func2(xkm1, xkm2, t) + func2(xkm1, xkm3, t) - func2(xkm2, xkm3, t);
        xkp = xkm1 - (2*calculate(xkm1, t))/(w + sqrt(w*w - 4*calculate(xkm1, t)*func3(xkm1, xkm2, xkm3, t)));
        xkm = xkm1 - (2*calculate(xkm1, t))/(w - sqrt(w*w - 4*calculate(xkm1, t)*func3(xkm1, xkm2, xkm3, t))); 
        xk=0;
        if(abs(xkm1-xkp)<abs(xkm1-xkm))
            xk = xkp;
        else
            xk = xkm;
        end
        xkm3=xkm2;
        xkm2 = xkm1;
        xkm1 = xk;
        iter = iter+1;
        fprintf('%5d \t %15d + %12d i \t %10d\n', iter, real(xk), imag(xk), calculate(xk,t)); 
    end
        fprintf('\nSolution =  %15d + %12d i \n\n', real(xk), imag(xk));
end

function f = func3(x, y, z, t)
    f = (func2(y, z, t)-func2(x, y, t))/(z-x);
end
function f = func2(x, y, t)
    f = (calculate(y, t)-calculate(x, t))/(y-x);
end

function f = calculate(x, t)
    f=0;
    if(t==1)
        f = x^3 - x- 2;
    end
    if(t==2)
        f = 1+2*x -tan(x);
    end
    if(t==3)
        f = x^2 + exp(x) -5;
    end
    if(t==4)
        f = x^2 - sin(x);

    end
    if(t==5)
        f = x^4 - 2*x^3 -2*1i*x^2 + 4*1i*x;
    end
    if(t==6)
        f = x - exp(x);
    end

end