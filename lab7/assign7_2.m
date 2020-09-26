clc;
clear all;
close all;
a=0;
b=1;
y0=1;

inc = [0.05 0.005];
ex = @(x) 2*exp(x)-x-1;
f = @(x, y)x+y;
sol = cell(1, 2);
for i=1:2
    sol{i} = solve_simple(a, b, inc(i), f, y0, ex);
end

max_error_table = create_max_table(2, sol);
disp('Using Simple Euler');
plot_and_print2(sol, 2, 1, 2, 1, max_error_table, 'Simple Euler: ');

sol2 = cell(1, 2);
for i=1:2
    sol2{i} = solve_modified(a, b, inc(i), f, y0, ex);
end

max_error_table2 = create_max_table(2, sol2);
fprintf('\n\n');
disp('Using Modified Euler');
plot_and_print2(sol2, 2, 1, 2, 4, max_error_table2, 'Modified Euler: ');

fprintf('\nEstimated value of y(1) using simple euler for h = 0.05: %d, absolute error : %d \n', sol{1}.Y(21), sol{1}.error(21));
fprintf('Estimated value of y(1) using simple euler for h = 0.005: %d, absolute error : %d \n', sol{2}.Y(201), sol{2}.error(201));
fprintf('\nEstimated value of y(1) using modified euler for h = 0.05: %d, absolute error : %d \n', sol2{1}.Y(21), sol2{1}.error(21));
fprintf('Estimated value of y(1) using modified euler for h = 0.005: %d, absolute error : %d \n', sol2{2}.Y(201), sol2{2}.error(201));
function plot_and_print2(sol, num, s1, s2, fig, max_error_table, w)
    figure(fig);
    for i=1:num
        subplot(s1, s2, i);
        hold on;
        plot(sol{i}.t, sol{i}.Y, 'color', 'y', 'LineWidth', 3);
        plot(sol{i}.t, sol{i}.exact,  'color', 'k', 'LineWidth', 1);
        title(strcat(w, 'Approximate Vs Actual for h = ',num2str(sol{i}.h)));
        legend('Approximate', 'Actual');
        hold off;
    end
    figure(fig+1);
    for i=1:num
        subplot(s1, s2, i);
        plot(sol{i}.t, sol{i}.error);
        title(strcat(w, 'Absolute error for h = ',num2str(sol{i}.h)));
    end
    for i=1:num
        table_summary(sol{i}.tab, sol{i}.h, sol{i}.n);
    end
    figure(fig+2);
    plot(log2(log2(max_error_table(:, 2))), max_error_table(:, 3));
    title(strcat(w, 'loglog(n) vs max error for all cases'));
    disp('Max error table for all cases');
    fprintf('\t%s\t\t\t\t\t%s\t\t\t\t\t%s\t\t\t\t%s\n', 'h', 'n', 'max error', 'log2(En/En+1)');
    for i=1:num
        fprintf('%10d\t\t%6d\t\t%25.5f\t\t%12.5f\n', max_error_table(i, :));
    end
end
function table_summary(tab, h, num)
    disp(strcat('table summary for h = ',num2str(h)));
    fprintf('\t\t%s\t\t\t\t\t%s\t\t\t\t%s\t\t\t\t\t%s\n', 't', 'Approximate Y', 'Exact Y', 'Error');
    for i=1:2
        fprintf('%13.5d\t\t%25.5f\t\t%10.5f\t\t%25.5f\n', tab(i, :));
    end
    for i=num-1:num
        fprintf('%13.5d\t\t%25.5f\t\t%10.5f\t\t%25.5f\n', tab(i, :));
    end
end


function max_error_table = create_max_table(num, sol)
    max_error_table = zeros(num, 4);
    for i=1:num
        max_error_table(i, 1) = sol{i}.h;
        max_error_table(i, 2) = sol{i}.n;
        max_error_table(i, 3) = sol{i}.max_error;
    end
    for i=1:num-1
        max_error_table(i, 4) = log2(max_error_table(i, 3)/max_error_table(i+1, 3));
    end
end

function sol = solve_simple(a, b, h, f, y0, ex)
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