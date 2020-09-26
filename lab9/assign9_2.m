clear all
clc
clf
close all

%problem (a)% 
h = 0.001;
a = -1;
b = 0;
alpha = -2;
beta = 1;
p = @(x)-x;
q = @(x)2-x+x;
r = @(x)-2-(2+x.*x).*exp(x);
ex = @(x)x.^2+x.*exp(x);
ans1 = solve_central(h, a, b, alpha, beta, p, q, r, ex);
print(ans1, 'problem (a) - ', 1);

%problem(b)%
h2 = 0.001;
a2 = 0;
b2 = 1;
alpha2 = -1;
beta2 = 2*sin(1);
p2 = @(x)(-x+x);
q2 = @(x)(-x);
r2 = @(x)(-3+x+x.*x-x.^3).*sin(x)-4.*x.*cos(x);
ex2 = @(x)(sin(pi*x)).^2;
ans2 = solve_central(h2, a2, b2, alpha2, beta2, p2, q2, r2, ex2);
print2(ans2, 'problem (b) - ', 3);

function X = solve_central(h, a, b, alpha, beta, p, q, r, ex)
    X.t = a:h:b;
    s = size(X.t, 2);

    u_prev = (-1/(h*h)-(p(X.t))/(2*h));
    u_cur = (2/(h*h)+q(X.t));
    u_next = (-1/(h*h)+(p(X.t))/(2*h));

    A = zeros(s, s);
    B = zeros(s, 1);
    for i=2:s-1
        A(i, i) = u_cur(i);
        B(i)= r(X.t(i));
    end
    A(1, 1)=-1/h;
    A(1, 2)= 1/h;
    A(s, s)=1/h;
    A(s, s-1)=-1/h;

    for i=2:s-1
        A(i, i-1) = u_prev(i);
        A(i, i+1) = u_next(i);
    end
    B(1)=alpha;
    B(s)=beta;
    X.u  = (A\B)';
    X.ex = ex(X.t);
    X.error = abs(X.ex-X.u);
    X.n=s;
end

function print(sol, str, f)
    disp(str);
    fprintf('\t%s\t%s\t%s\t%s\n', 'x', 'Approximate y(x)', 'Exact y(x)', 'Error');
    for i=1:100:sol.n
        fprintf('%15e\t%15e\t%15e\t%15e\n', sol.t(i), sol.u(i), sol.ex(i), sol.error(i));
    end
    figure(f);
    hold on;
    plot(sol.t, sol.u, 'Linewidth', 3, 'color', 'y');
    plot(sol.t, sol.ex, 'Linewidth', 1, 'color', 'k');
    legend('Approximate', 'Exact');
    title(strcat(str, 'Approximate Vs Exact Solution'));
    hold off;
    
    figure(f+1);
    plot(sol.t, sol.error);
    title(strcat(str, 'Error in the solution y(x)'));
end
function print2(sol, str, f)
    disp(str);
    fprintf('\t%s\t%s\n', 'x', 'Approximate y(x)');
    for i=1:100:sol.n
        fprintf('%15e\t%15e\n', sol.t(i), sol.u(i));
    end
    figure(f);
    hold on;
    plot(sol.t, sol.u, 'Linewidth', 2, 'color', 'm');
    title(strcat(str, 'Approximate y(x)'));
    hold off;
end