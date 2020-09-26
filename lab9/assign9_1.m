clear all
clc
clf
close all
%problem (a)% 
h = 0.001;
a = 0;
b = 1;
alpha = 1;
beta = 3;
p = @(x)-x-1;
q = @(x)cos(x);
r = @(x)-exp(x);

sol = solve(h, a, b, alpha, beta, p, q, r);
print(sol, 'problem (a) - ', 1);
%problem(b)%

h2 = 0.001;
a2 = 0;
b2 = 1;
alpha2 = 0;
beta2 = 0;
p2 = @(x)(-2-x+x);
q2 = @(x)(-1-x+x);
r2 = @(x)-x;

sol2 = solve(h, a2, b2, alpha2, beta2, p2, q2, r2);
print(sol2, 'problem (b) - ', 2);

function sol = solve(h, a, b, alpha, beta, p, q, r)
    sol.t = a:h:b;
    sol.n = size(sol.t, 2);
    sol.central = solve_central(sol.t, sol.n, h, alpha, beta, p, q, r);
    sol.forward = solve_forward(sol.t, sol.n, h, alpha, beta, p, q, r);
    sol.backward = solve_backward(sol.t, sol.n, h, alpha, beta, p, q, r);
    
end


function X = solve_central(t, s, h, alpha, beta, p, q, r)
    u_prev = (-1/(h*h)-(p(t))/(2*h));
    u_cur = (2/(h*h)+q(t));
    u_next = (-1/(h*h)+(p(t))/(2*h));

    A = zeros(s, s);
    B = zeros(s, 1);
    for i=2:s-1
        A(i, i) = u_cur(i);
        B(i)= r(t(i));
    end
    A(1, 1)=1;
    A(s, s)=1;

    for i=2:s-1
        A(i, i-1) = u_prev(i);
        A(i, i+1) = u_next(i);
    end
    B(1)=alpha;
    B(s)=beta;
    X  = A\B;

end

function X = solve_forward(t, s, h, alpha, beta, p, q, r)
    u_prev = (-1/(h*h));
    u_cur = (2/(h*h)+q(t)-(p(t))/(h));
    u_next = (-1/(h*h)+(p(t))/(h));

    A = zeros(s, s);
    B = zeros(s, 1);
    for i=2:s-1
        A(i, i) = u_cur(i);
        B(i)= r(t(i));
    end
    A(1, 1)=1;
    A(s, s)=1;

    for i=2:s-1
        A(i, i-1) = u_prev;
        A(i, i+1) = u_next(i);
    end
    B(1)=alpha;
    B(s)=beta;
    X  = A\B;

end

function X = solve_backward(t, s, h, alpha, beta, p, q, r)


    u_prev = (-1/(h*h)-(p(t))/(h));
    u_cur = (2/(h*h)+q(t)+(p(t))/(h));
    u_next = (-1/(h*h));

    A = zeros(s, s);
    B = zeros(s, 1);
    for i=2:s-1
        A(i, i) = u_cur(i);
        B(i)= r(t(i));
    end
    A(1, 1)=1;
    A(s, s)=1;

    for i=2:s-1
        A(i, i-1) = u_prev(i);
        A(i, i+1) = u_next;
    end
    B(1)=alpha;
    B(s)=beta;
    X  = A\B;

end

function print(sol, str, f)
    disp(str);
    fprintf('\t%s\t%s\t%s\t%s\n', 'x', 'y(x): Backward Case', 'y(x): Forward Case', 'y(x): Central Case');
    for i=1:100:sol.n
        fprintf('%12e\t%15e\t\t%15e\t\t%15e\n', sol.t(i), sol.backward(i), sol.forward(i), sol.central(i));
    end
    figure(f);
    hold on;
    
    plot(sol.t, sol.backward, 'Linewidth', 3, 'color', 'y');
    plot(sol.t, sol.forward, 'Linewidth', 4, 'color', 'm', 'Linestyle', '--');
    plot(sol.t, sol.central, 'Linewidth', 1, 'color', 'k');
    legend('Using Backward Difference', '(Using Forward Difference)', '(Using Central Difference)');
    title(strcat(str, 'Approximate y(x) using different finite difference schemes for first order derivative'));
    hold off;
    fprintf('\n\n');
end