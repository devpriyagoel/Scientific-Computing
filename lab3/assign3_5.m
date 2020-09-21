clear all
clc
clf
close all
syms x
func = x^5 - 2*x^4 - 2*x^3+8*x^2-7*x+2;
gunc = x^5-8*x^4+25*x^3-38*x*x+28*x-8;
f = @(x) x^5 - 2*x^4 - 2*x^3+8*x^2-7*x+2;
g = @(x) x^5-8*x^4+25*x^3-38*x*x+28*x-8;
dfunc = diff(func)
dgunc = diff(gunc)
df = @(x) 5*x^4 - 8*x^3 - 6*x^2 + 16*x - 7;
dg = @(x) 5*x^4 - 32*x^3 + 75*x^2 - 76*x + 28;
t = 1.3;
p = 4;
modified_newton(f, df, p, t, 1)
t = 1.3;
p = 2;
modified_newton(g, dg, p, t, 1)
t = 3;
p = 3;
modified_newton(g, dg, p, t, 2)

function modified_newton(func, df, p, t, root)
    del = 1;
    iter=0;
    e=0;
    eprev=0;
    fprintf('iter\t\tx \t\t f(x)\n');
    while(abs(t-root)>1e-5)
        j = df(t);
        f = func(t);
        eprev = e;
        e = abs(root-t);
        fprintf('%d\t%13d\t%13d\t%13d\t%13d\t%13d\n', iter, t, f, j, e/eprev, log(e/eprev));
        del = f/j;
        t = t-p*del;
        iter = iter+1;
        
    end
    eprev = e;
    e = abs(root-t);
    fprintf('%d\t%13d\t%13d\t%13d\t%13d\t%13d\n', iter, t, func(t), df(t), e/eprev, log(e/eprev));
    fprintf('\n Solution : %14d\n', t);
end

