clc;
clear all;
close all;
a=0;
b=10;
y0=1;

inc = [0.01 0.001 0.005 0.0025];
ex = @(x) sin(x)+cos(x);
f = @(x, y)-5*y + 6*cos(x)+4*sin(x);
sol = cell(1, 6);
for i=1:4
    sol{i} = solve_explicit(a, b, inc(i), f, y0, ex, -5);
end

f2 = @(x, y)5*y-4*cos(x)-6*sin(x);
sol{5} = solve_explicit(a, b, inc(1), f2, y0, ex, 5);
sol{6} = solve_explicit(a, b, inc(2), f2, y0, ex, 5);
max_error_table = create_max_table(6, sol);
% plot_and_print_with_uitable(sol, 6, 1, max_error_table);
disp('Using Explicit Euler');
plot_and_print2(sol, 6, 3, 2, 1, max_error_table, 'Explicit Euler: ');

sol2 = cell(1, 6);
f = @(x, y)-5*y + 6*cos(x)+4*sin(x);
for i=1:4
    sol2{i} = solve_implicit(a, b, inc(i), f, y0, ex, -5);
end

f2 = @(x, y)5*y-4*cos(x)-6*sin(x);
sol2{5} = solve_implicit(a, b, inc(1), f2, y0, ex, 5);
sol2{6} = solve_implicit(a, b, inc(2), f2, y0, ex, 5);
max_error_table2 = create_max_table(6, sol2);
fprintf('\n\n');
% plot_and_print_with_uitable(sol2, 6, 1, max_error_table2);
disp('Using Implicit Euler');
plot_and_print2(sol2, 6, 3, 2, 4, max_error_table2, 'Implicit Euler: ');


function plot_and_print2(sol, num, s1, s2, fig, max_error_table, w)
    figure(fig);
    for i=1:num
        subplot(s1, s2, i);
        hold on;
        plot(sol{i}.t, sol{i}.Y, 'color', 'y', 'LineWidth', 3);
        plot(sol{i}.t, sol{i}.exact,  'color', 'k', 'LineWidth', 1);
        title(strcat(w, 'Approximate Vs Actual for lambda = ', int2str(sol{i}.lambda),', h = ',num2str(sol{i}.h)));
        legend('Approximate', 'Actual');
        hold off;
    end
    figure(fig+1);
    for i=1:num
        subplot(s1, s2, i);
        plot(sol{i}.t, sol{i}.error);
        title(strcat(w, 'Absolute error for lambda = ', int2str(sol{i}.lambda),', h = ',num2str(sol{i}.h)));
    end
    for i=1:num
        table_summary(sol{i}.tab, sol{i}.lambda, sol{i}.h, sol{i}.n);
    end
    figure(fig+2);
    plot(log2(log2(max_error_table(:, 3))), max_error_table(:, 4));
    title(strcat(w, 'loglog(n) vs max error for all 6 cases'));
    disp('Max error table for all 6 cases');
    fprintf('%s\t\t%s\t\t\t\t\t%s\t\t\t\t\t%s\t\t\t\t%s\n', 'lambda', 'h', 'n', 'max error', 'log2(En/En+1)');
    for i=1:num
        fprintf('%3d\t\t%10d\t\t%6d\t\t%25.5f\t\t%12.5f\n', max_error_table(i, :));
    end
end
function table_summary(tab, lambda, h, num)
    disp(strcat('table summary for lambda = ', int2str(lambda),', h = ',num2str(h)));
    fprintf('\t\t%s\t\t\t\t\t%s\t\t\t\t%s\t\t\t\t\t%s\n', 't', 'Approximate Y', 'Exact Y', 'Error');
    for i=1:2
        fprintf('%13.5d\t\t%25.5f\t\t%10.5f\t\t%25.5f\n', tab(i, :));
    end
    for i=num-1:num
        fprintf('%13.5d\t\t%25.5f\t\t%10.5f\t\t%25.5f\n', tab(i, :));
    end
end


function max_error_table = create_max_table(num, sol)
    max_error_table = zeros(num, 5);
    for i=1:num
        max_error_table(i, 1) = sol{i}.lambda;
        max_error_table(i, 2) = sol{i}.h;
        max_error_table(i, 3) = sol{i}.n;
        max_error_table(i, 4) = sol{i}.max_error;
    end
    for i=1:num-1
        max_error_table(i, 5) = log2(max_error_table(i, 4)/max_error_table(i+1, 4));
    end
end
% function plot_and_print_with_uitable(sol, num, fig, max_error_table)
%     VarNames = {'t', 'Approximate Y', 'Exact Y', 'Error'};
%     for i=1:num
%         uitable(figure(fig+i-1),'data',sol{i}.tab,'columnname',VarNames);
%         subplot(2, 2, 2);
%         hold on;
%         plot(sol{i}.t, sol{i}.Y, 'color', 'y', 'LineWidth', 3);
%         plot(sol{i}.t, sol{i}.exact,  'color', 'k', 'LineWidth', 1);
%         title(strcat('Approximate Vs Actual for lambda = ', int2str(sol{i}.lambda),', h = ',num2str(sol{i}.h)));
%         legend('Approximate', 'Actual');
%         hold off;
%         subplot(2, 2, 4);
%         plot(sol{i}.t, sol{i}.error);
%         title(strcat('Absolute error for lambda = ', int2str(sol{i}.lambda),', h = ',num2str(sol{i}.h)));
%     end
%     VarNames = ['lambda', 'h', 'n', 'max error', 'log2(En/En+1)'];
%     uitable(figure(fig+num),'data',max_error_table,'columnname',VarNames);
%     subplot(1, 2, 2);
%     plot(log2(log2(max_error_table(:, 3))), max_error_table(:, 4));
%     title('loglog(max_error) for all 6 cases');
% end
function sol = solve_explicit(a, b, h, f, y0, ex, lambda)
    sol.t = a:h:b;
    sol.h = h;
    sol.lambda = lambda;
    sol.n = size(sol.t, 2);
    sol.exact = ex(sol.t);
    sol = explicit(sol, h, f, y0);
    sol.tab = [transpose(sol.t),transpose(sol.Y),transpose(sol.exact),transpose(sol.error)];
end
function sol = solve_implicit(a, b, h, f, y0, ex, lambda)
    sol.t = a:h:b;
    sol.h = h;
    sol.lambda = lambda;
    sol.n = size(sol.t, 2);
    sol.exact = ex(sol.t);
    sol = implicit(sol, h, f, y0);
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