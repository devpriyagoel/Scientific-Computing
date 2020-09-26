clc;
clear all;
close all;

q.a = 0;
q.b = 1;
q.h = 2^(-7);
q.y0 = 0;
q.ex = @(t)t.*asin(t);
q.f = @(t, y) y./t+ t.*sec(y./t);
q.n = (q.b-q.a)./q.h +1;
q.num = 1;

sol = Runge_Kutta_fourth_ord(q);

function s = Runge_Kutta_fourth_ord(q)
    s = cell(1, q.num);
    for i=1:q.num
        s{i} = zeros(q.n(i), 4);
        s{i}(:, 1) = q.a:q.h(i):q.b;
        s{i}(:, 3) = q.ex(s{i}(:, 1));
        s{i}(1, 2) = q.y0;
        for j=2:2   %special case: defined : f(0, 0) = 0%
            f1 = 0;
            f2 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f1/2);
            f3 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f2/2);
            f4 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i), s{i}(j-1, 2)+f3);
            s{i}(j, 2) = s{i}(j-1, 2)+ (f1+2*f2+2*f3+f4)/6;
        end
        for j=3:q.n(i)
            f1 = q.h(i)*q.f(s{i}(j-1, 1), s{i}(j-1, 2));
            f2 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f1/2);
            f3 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f2/2);
            f4 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i), s{i}(j-1, 2)+f3);
            s{i}(j, 2) = s{i}(j-1, 2)+ (f1+2*f2+2*f3+f4)/6;
        end
        s{i}(:, 4) = abs(s{i}(:, 2)-s{i}(:, 3));
        figure(i);
        subplot(1, 2, 1);
        hold on;
        plot(s{i}(:, 1), s{i}(:, 2),'color', 'y', 'LineWidth', 3);
        plot(s{i}(:, 1), s{i}(:, 3),'color', 'k', 'LineWidth', 1);
        title(strcat('4th order Runge Kutta: Approximate Vs Actual for h=', num2str(q.h(i))));
        legend('Approximate', 'Actual');
        hold off;
        subplot(1, 2, 2);
        plot(s{i}(:, 1), s{i}(:, 4),'LineWidth', 2);
        title(strcat('4th order Runge Kutta: Absolute Error for h=', num2str(q.h(i))));
        disp(strcat('4th order Runge Kutta: table for h = ',num2str(q.h(i))));
        fprintf('\t\t%s\t\t\t%s\t\t\t%s\t\t\t%s\n', 't', 'Approximate Y', 'Exact Y', 'Error');
        for j=1:q.n(i)
            fprintf('%12d\t\t%12.5f\t\t%10.5f\t\t%10.5f\n', s{i}(j, :));
        end
    end
end