clear all
clc
clf
close all
syms x
f = @(x) cos(x+sqrt(2))+x*(x/2+sqrt(2));
g = @(x) exp(6*x)+3*(log(2))^2*exp(2*x)-log(8)*exp(4*x)-(log(2))^3;
syms x
Jf = jacobian([f],[x])
Jg = jacobian([g],[x])
diff = @(x) x - sin(x + 2^(1/2)) + 2^(1/2);
difg = @(x) (3245652021676565*exp(2*x))/1125899906842624 - (4682486076124019*exp(4*x))/562949953421312 + 6*exp(6*x);

t = -1.5;
t = newton(f, diff, t);
p = multiplicity(Jf, t)
t = -1.5;
modified_newton(f, diff, p, t)
t = -0.5;
t = newton(g, difg, t);
p = multiplicity(Jg, t)
t = -0.5;
modified_newton(g, difg, p, t)
function root = newton(func, diff, t)
    del = 1;
    iter=0;
    fprintf('iter\t\tx \t\t f(x)\n');
    while(abs(del)>1e-6)
        j = diff(t);
        f = func(t);
        fprintf('%d\t%d\t%d\n', iter, t, f);
        del = f/j;
        t = t-del;
        iter = iter+1;
    end
    fprintf('%d\t%d\t%d\n', iter, t, func(t));
    fprintf('\n Solution : %14d\n', t);
    root = t;
end
function modified_newton(func, diff, p, t)
    del = 1;
    iter=0;
    fprintf('iter\t\tx \t\t f(x)\n');
    while(abs(del)>1e-6)
        j = diff(t);
        f = func(t);
        fprintf('%d\t%d\t%d\n', iter, t, f);
        del = f/j;
        t = t-p*del;
        iter = iter+1;
        
    end
    fprintf('%d\t%d\t%d\n', iter, t, func(t));
    fprintf('\n Solution : %14d\n', t);
end
function p = multiplicity(F,root)
    df = matlabFunction(F);
    F = diff(F);
    p = 1;
    while(abs(df(root))<1e-3)
        F = diff(F);
        df = matlabFunction(F);
        p = p+1;
    end
end
