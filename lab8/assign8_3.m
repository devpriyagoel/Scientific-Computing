clc;
clear all;
close all;

q.a = 0;
q.b = 0.2;
q.h = 1e-5;
q.y0 = 0;
q.f = @(t, y) (6.22*(10^(-19))).*((2000-y/2)^4).*((3000-3*y/4)^3);
q.n = (q.b-q.a)./q.h +1;
q.num = 1;

sol = Runge_Kutta_fourth_ord(q);

function s = Runge_Kutta_fourth_ord(q)
    s = cell(1, q.num);
    for i=1:q.num
        s{i} = zeros(q.n(i), 2);
        s{i}(:, 1) = q.a:q.h(i):q.b;
        s{i}(1, 2) = q.y0;
        for j=2:q.n(i)
            f1 = q.h(i)*q.f(s{i}(j-1, 1), s{i}(j-1, 2));
            f2 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f1/2);
            f3 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f2/2);
            f4 = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i), s{i}(j-1, 2)+f3);
            s{i}(j, 2) = s{i}(j-1, 2)+ (f1+2*f2+2*f3+f4)/6;
        end
        figure(i);
        plot(s{i}(:, 1), s{i}(:, 2),'color', 'k', 'LineWidth', 2);
        title('Approximate Amount of KOH formed vs time');
        disp('Approximate amount of KOH formed at 0.2 seconds:');
        disp(s{i}(q.n(i), 2));
    end
end
