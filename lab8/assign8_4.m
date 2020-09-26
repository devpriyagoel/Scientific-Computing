clc;
clear all;
close all;

q.a = 1;
q.b = 3;
q.h = 0.002;
q.y0 = 0;
q.ex = @(t)t.*tan(log(t));
q.f = @(t, y) 1+y./t+ (y./t).^2;
q.n = (q.b-q.a)./q.h +1;
q.num = 1;

r.a = 0;
r.b = 1;
r.h = 0.005;
r.y0 = 1;
r.ex = @(t)(sin(t)-cos(t)+3*exp(-t))/2;
r.f = @(t, y) sin(t)-y;
r.n = (r.b-r.a)./r.h +1;
r.num = 1;

sol1 = Adams_Bashforth(q, 0, 1);
sol2 = Adams_Bashforth(q, 1, 2);
sol3 = Adams_Moulton(q, 0, 3);
sol4 = Adams_Moulton(q, 1, 4);
sol5 = Adams_Bashforth(r, 0, 5);
sol6 = Adams_Bashforth(r, 1, 6);
sol7 = Adams_Moulton(r, 0, 7);
sol8 = Adams_Moulton(r, 1, 8);

function s = Adams_Bashforth(q, flag, fig)
    s = cell(1, q.num);
    for i=1:q.num
        s{i} = zeros(q.n(i), 4);
        s{i}(:, 1) = q.a:q.h(i):q.b;
        s{i}(:, 3) = q.ex(s{i}(:, 1));
        s{i}(1, 2) = q.y0;
        f = zeros(1, 4);
        if(flag==1)
            for j=2:4
                f(1) = q.h(i)*q.f(s{i}(j-1, 1), s{i}(j-1, 2));
                f(2) = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f(1)/2);
                f(3) = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f(2)/2);
                f(4) = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i), s{i}(j-1, 2)+f(3));
                s{i}(j, 2) = s{i}(j-1, 2)+ (f(1)+2*f(2)+2*f(3)+f(4))/6;
            end
        else
            for j=2:4
                s{i}(j, 2) = q.ex(s{i}(j, 1));
            end
        end
        
        for j=1:4
            f(j)=q.f(s{i}(j, 1), s{i}(j, 2));
        end
        for j=5:q.n(i)
            s{i}(j, 2) = s{i}(j-1, 2)+(q.h(i)/24)*(55*f(4)-59*f(3)+37*f(2)-9*f(1));
            for k=1:3
                f(k)=f(k+1);
            end
            f(4)=q.f(s{i}(j, 1), s{i}(j, 2));
        end
        s{i}(:, 4) = abs(s{i}(:, 2)-s{i}(:, 3));
        figure(i+fig-1);
        subplot(1, 2, 1);
        hold on;
        plot(s{i}(:, 1), s{i}(:, 2),'color', 'y', 'LineWidth', 3);
        plot(s{i}(:, 1), s{i}(:, 3),'color', 'k', 'LineWidth', 1);
        title(strcat('Adams Bashforth: Approximate Vs Actual for h=', num2str(q.h(i))));
        legend('Approximate', 'Actual');
        hold off;
        subplot(1, 2, 2);
        plot(s{i}(:, 1), s{i}(:, 4),'LineWidth', 2);
        title(strcat('Adams Bashforth: Absolute Error for h=', num2str(q.h(i))));
        disp(strcat('Adams Bashforth: table for h = ',num2str(q.h(i))));
        fprintf('\t\t%s\t\t\t%s\t\t\t%s\t\t\t%s\n', 't', 'Approximate Y', 'Exact Y', 'Error');
        m = 50;
        for j=1:m:q.n(i)
            fprintf('%12d\t\t%12.5e\t\t%10.5f\t\t%11.9e\n', s{i}(j, :));
        end
    end
end
function s = Adams_Moulton(q, flag, fig)
    s = cell(1, q.num);
    tol = 1e-15;
    for i=1:q.num
        s{i} = zeros(q.n(i), 4);
        s{i}(:, 1) = q.a:q.h(i):q.b;
        s{i}(:, 3) = q.ex(s{i}(:, 1));
        s{i}(1, 2) = q.y0;
        f = zeros(1, 4);
        if(flag==1)
            for j=2:4
                f(1) = q.h(i)*q.f(s{i}(j-1, 1), s{i}(j-1, 2));
                f(2) = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f(1)/2);
                f(3) = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i)/2, s{i}(j-1, 2)+f(2)/2);
                f(4) = q.h(i)*q.f(s{i}(j-1, 1)+q.h(i), s{i}(j-1, 2)+f(3));
                s{i}(j, 2) = s{i}(j-1, 2)+ (f(1)+2*f(2)+2*f(3)+f(4))/6;
            end
        else
            for j=2:4
                s{i}(j, 2) = q.ex(s{i}(j, 1));
            end
        end
        
        for j=1:4
            f(j)=q.f(s{i}(j, 1), s{i}(j, 2));
        end
        for j=5:q.n(i)
            s{i}(j, 2) = s{i}(j-1, 2)+(q.h(i)/24)*(55*f(4)-59*f(3)+37*f(2)-9*f(1));
            for k=1:3
                f(k)=f(k+1);
            end
            f(4)=q.f(s{i}(j, 1), s{i}(j, 2));
            prev = s{i}(j, 2)-1;
            l=1;
            while(abs(s{i}(j, 2)-prev)<tol&&l<30)
                prev = s{i}(j, 2);
                s{i}(j, 2) = s{i}(j-1, 2)+(q.h(i)/24)*(9*f(4)+19*f(3)-5*f(2)+f(1));
                f(4)=q.f(s{i}(j, 1), s{i}(j, 2));
                l = l+1;
            end
        end
        s{i}(:, 4) = abs(s{i}(:, 2)-s{i}(:, 3));
        figure(i+fig-1);
        subplot(1, 2, 1);
        hold on;
        plot(s{i}(:, 1), s{i}(:, 2),'color', 'y', 'LineWidth', 3);
        plot(s{i}(:, 1), s{i}(:, 3),'color', 'k', 'LineWidth', 1);
        title(strcat('Adams Moulton: Approximate Vs Actual for h=', num2str(q.h(i))));
        legend('Approximate', 'Actual');
        hold off;
        subplot(1, 2, 2);
        plot(s{i}(:, 1), s{i}(:, 4),'LineWidth', 2);
        title(strcat('Adams Moulton: Absolute Error for h=', num2str(q.h(i))));
        disp(strcat('Adams Moulton: table for h = ',num2str(q.h(i))));
        fprintf('\t\t%s\t\t\t%s\t\t\t%s\t\t\t%s\n', 't', 'Approximate Y', 'Exact Y', 'Error');
        m = 50;
        for j=1:m:q.n(i)
            fprintf('%12d\t\t%12.5e\t\t%10.5f\t\t%11.9e\n', s{i}(j, :));
        end
    end
end
