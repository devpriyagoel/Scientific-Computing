clear all
clc
clf
close all
%problem (a)% 
h = 0.001;
a = 0;
b = 1;
alpha = 1;
beta = 1;
coeff = [(1-1/h) 1/h 0 1];
p = @(x)-x+x;
q = @(x)x;
r = @(x)-1+x-x;

sol = solve(h, a, b, alpha, beta, p, q, r, coeff);
print(sol, 'problem (a) - ', 1);
%problem(b)%

h2 = 0.001;
a2 = 0;
b2 = 1;
alpha2 = 1;
beta2 = 1;
coeff2 = [(-1-1/h) 1/h -1/h (1+1/h)];
p2 = @(x)(3+x-x);
q2 = @(x)(-2+x-x);
r2 = @(x)-2+x-x;

sol2 = solve(h, a2, b2, alpha2, beta2, p2, q2, r2, coeff2);
print(sol2, 'problem (b) - ', 2);

function sol = solve(h, a, b, alpha, beta, p, q, r, coeff)
    sol.t = a:h:b;
    sol.n = size(sol.t, 2);
    sol.central = solve_central(sol.t, sol.n, h, alpha, beta, p, q, r, coeff);
    
end


function X = solve_central(t, s, h, alpha, beta, p, q, r, coeff)
    u_prev = (-1/(h*h)-(p(t))/(2*h));
    u_cur = (2/(h*h)+q(t));
    u_next = (-1/(h*h)+(p(t))/(2*h));

    A = zeros(s, s);
    B = zeros(s, 1);
    for i=2:s-1
        A(i, i) = u_cur(i);
        B(i)= r(t(i));
    end
    A(1, 1)=coeff(1);
    A(1, 2)=coeff(2);
    A(s, s-1)=coeff(3);
    A(s, s)=coeff(4);
    for i=2:s-1
        A(i, i-1) = u_prev(i);
        A(i, i+1) = u_next(i);
    end
    B(1)=alpha;
    B(s)=beta;
    X  = A\B;

end


function print(sol, str, f)
    disp(str);
    fprintf('\t%s\t%s\n', 'x', 'y(x): Approximate');
    for i=1:100:sol.n
        fprintf('%12e\t%15e\n', sol.t(i), sol.central(i));
    end
    figure(f);
    plot(sol.t, sol.central, 'Linewidth', 2, 'color', 'k');
    title(strcat(str, 'Approximate y(x) for mixed conditions'));
    fprintf('\n\n');
end