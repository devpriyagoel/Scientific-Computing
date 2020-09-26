clc;
clear all;
close all;
a=[0 1 1 0];
b=[1 3 3 1];
y0=[1 -2 0 1];

inc = [0.01 0.02 0.002 0.005];
ex = {@(x)(2.*x + 1)./(x.*x + 1), @(x)(2.*x)./(1-2.*x), @(x)x.*tan(log(x)), @(x)log(exp(x) + exp(1) -1)};  
f = {@(x, y)(2 - 2*x.*y)./(1 + x.*x), @(x, y)(y^2 + y)/x , @(x, y)1 + (y/x) + (y/x)^2, @(x, y)exp(x-y)};
sol = cell(4, 4);
for i=1:4
    sol{1, i} = solve_explicit(a(i), b(i), inc(i), f{i}, y0(i), ex{i});
    sol{2, i} = solve_implicit(a(i), b(i), inc(i), f{i}, y0(i), ex{i});
    sol{3, i} = solve_modified(a(i), b(i), inc(i), f{i}, y0(i), ex{i});
    sol{4, i} = solve_trapezoid(a(i), b(i), inc(i), f{i}, y0(i), ex{i});
end

plot_and_print(sol);
function plot_and_print(sol)
    head = {'problem (a)', 'problem (b)', 'problem (c)', 'problem (d)'};
    head2 = {'Explicit Euler: ', 'Implicit Euler: ', 'Modified Euler: ', 'Trapezoidal Euler: '};
    for i=1:4
        figure(i);
        for j=1:4
            subplot(2, 4, j);
            hold on;
            plot(sol{j, i}.t, sol{j, i}.Y, 'color', 'y', 'LineWidth', 3);
            plot(sol{j, i}.t, sol{j, i}.exact,  'color', 'k', 'LineWidth', 1);
            title(strcat(head2{j}, ' ',head{i}));
            legend('Approximate', 'Actual');
            hold off;
        end
        for j=1:4
            subplot(2, 4, j+4);
            plot(sol{j, i}.t, sol{j, i}.error);
            title('Absolute error');
        end
        for j=1:4
            table_summary(sol{j, i}.tab, head{i}, head2{j}, sol{j, i}.n);
        end
    end
    
end
function table_summary(tab, s1, s2, num)
    disp(strcat(s2, ' table summary for : ', s1));
    fprintf('\t\t%s\t\t\t%s\t\t%s\t\t%s\n', 't', 'Approximate Y', 'Exact Y', 'Error');
    for i=1:2
        fprintf('%13.5d\t\t%10.5f\t\t%10.5f\t\t%10.5f\n', tab(i, :));
    end
    for i=num-1:num
        fprintf('%13.5d\t\t%10.5f\t\t%10.5f\t\t%10.5f\n', tab(i, :));
    end
end

function sol = solve_explicit(a, b, h, f, y0, ex)
    sol.t = a:h:b;
    sol.h = h;
    sol.n = size(sol.t, 2);
    sol.exact = ex(sol.t);
    sol = explicit(sol, h, f, y0);
    sol.tab = [transpose(sol.t),transpose(sol.Y),transpose(sol.exact),transpose(sol.error)];
end
function sol = solve_modified(a, b, h, f, y0, ex)
    sol.t = a:h:b;
    sol.h = h;
    sol.n = size(sol.t, 2);
    sol.exact = ex(sol.t);
    sol = modified(sol, h, f, y0);
    sol.tab = [transpose(sol.t),transpose(sol.Y),transpose(sol.exact),transpose(sol.error)];
end

function sol = explicit(sol, h, f, y0)
    sol.Y = zeros(1, sol.n);
    sol.Y(1) = y0;
    for i=1:sol.n-1
        sol.Y(i+1) = sol.Y(i)+h*f(sol.t(i), sol.Y(i));
    end
    sol.error = abs(sol.Y-sol.exact);
    sol.max_error = max(sol.error); 
end
function sol = modified(sol, h, f, y0)
    sol.Y = zeros(1, sol.n);
    sol.Y(1) = y0;
    for i=1:sol.n-1
        sol.Y(i+1) = sol.Y(i)+(h/2)*(f(sol.t(i), sol.Y(i))+f(sol.t(i+1), sol.Y(i)+h*(f(sol.t(i), sol.Y(i)))));
    end
    sol.error = abs(sol.Y-sol.exact);
    sol.max_error = max(sol.error); 
end
function sol = solve_trapezoid(a, b, h, f, y0, ex)
    sol.t = a:h:b;
    sol.h = h;
    sol.n = size(sol.t, 2);
    sol.exact = ex(sol.t);
    sol = trapezoid(sol, h, f, y0);
    sol.tab = [transpose(sol.t),transpose(sol.Y),transpose(sol.exact),transpose(sol.error)];
end
function sol = trapezoid(sol, h, f, y0)
    sol.Y = zeros(1, sol.n);
    sol.Y(1) = y0;
    for i=1:sol.n-1
        it1=sol.Y(i);
        it2=1+it1;
        while(abs(it2-it1)>1e-4)
            it2 = sol.Y(i)+(h/2)*(f(sol.t(i), sol.Y(i))+f(sol.t(i+1), it1));
            it1 = it2;
        end
        sol.Y(i+1) = it2;
    end
    sol.error = abs(sol.Y-sol.exact);
    sol.max_error = max(sol.error); 
end
function sol = solve_implicit(a, b, h, f, y0, ex)
    sol.t = a:h:b;
    sol.h = h;
    sol.n = size(sol.t, 2);
    sol.exact = ex(sol.t);
    sol = implicit(sol, h, f, y0);
    sol.tab = [transpose(sol.t),transpose(sol.Y),transpose(sol.exact),transpose(sol.error)];
end
function sol = implicit(sol, h, f, y0)
    sol.Y = zeros(1, sol.n);
    sol.Y(1) = y0;
    for i=1:sol.n-1
        it1=sol.Y(i);
        it2=1+it1;
        while(abs(it2-it1)>1e-4)
            it2 = sol.Y(i)+h*f(sol.t(i), it1);
            it1 = it2;
        end
        sol.Y(i+1) = it2;
    end
    sol.error = abs(sol.Y-sol.exact);
    sol.max_error = max(sol.error); 
end