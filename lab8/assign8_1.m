clc;
clear all;
close all;

q.a = 0;
q.b = 5;
q.h = [0.05 0.025 0.0125 0.00625];
q.y0 = 0;
q.ex = @(t)t.^1.1;
q.f = @(t, y) -y+ (t.^0.1).*(1.1+t);
q.n = (q.b-q.a)./q.h +1;
q.num = 4;

sol = Runge_Kutta_sec_ord(q);

max_error_table = zeros(q.num, 3);
max_error_table(:, 1) = q.n;
for i=1:q.num
    max_error_table(i, 2) = max(sol{i}(:, 4));
end
for i=1:q.num-1
    max_error_table(i, 3) = log2(max_error_table(i, 2)/max_error_table(i+1, 2));
end
figure(q.num+1);
plot(log2(log2(max_error_table(:, 1))), max_error_table(:, 2));
title('2nd order Runge Kutta: loglog(n) vs max error');
disp('2nd order Runge Kutta: Max error table');
fprintf('\t%s\t\t\t\t%s\t\t%s\n', 'n', 'max error', 'log2(En/E2n)');
for i=1:q.num
    fprintf('%6d\t\t%15.5f\t\t%12.5f\n', max_error_table(i, :));
end


function s = Runge_Kutta_sec_ord(q)
    s = cell(1, q.num);
    for i=1:q.num
        s{i} = zeros(q.n(i), 4);
        s{i}(:, 1) = q.a:q.h(i):q.b;
        s{i}(:, 3) = q.ex(s{i}(:, 1));
        s{i}(1, 2) = q.y0;
        for j=2:q.n(i)
            f1 = q.h(i)*q.f(s{i}(j-1, 1), s{i}(j-1, 2));
            f2 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f1/2);
            s{i}(j, 2) = s{i}(j-1, 2)+f2;
        end
        s{i}(:, 4) = abs(s{i}(:, 2)-s{i}(:, 3));
        figure(i);
        subplot(1, 2, 1);
        hold on;
        plot(s{i}(:, 1), s{i}(:, 2),'color', 'y', 'LineWidth', 3);
        plot(s{i}(:, 1), s{i}(:, 3),'color', 'k', 'LineWidth', 1);
        title(strcat('2nd order Runge Kutta: Approximate Vs Actual for h=', num2str(q.h(i))));
        legend('Approximate', 'Actual');
        hold off;
        subplot(1, 2, 2);
        plot(s{i}(:, 1), s{i}(:, 4),'LineWidth', 2);
        title(strcat('2nd order Runge Kutta: Absolute Error for h=', num2str(q.h(i))));
        disp(strcat('2nd order Runge Kutta: table for h = ',num2str(q.h(i))));
        fprintf('\t%s\t\t\t%s\t\t%s\t\t\t%s\n', 't', 'Approximate Y', 'Exact Y', 'Error');
        m = 1/q.h(i);
        for j=m+1:m:q.n(i)
            fprintf('%5d\t\t%12.5f\t\t%10.5f\t\t%11.5f\n', s{i}(j, :));
        end
    end
end