clear all
clc
clf
close all

obj.x = 1:0.5:6;
obj.f = @(x) atan(x);
obj.y = obj.f(obj.x);
obj.n = size(obj.x, 2);

query.v = 0:0.25:8;
query.fd = newtons_forward_difference(obj, query.v);
query.exact = obj.f(query.v);
query.diff = abs(query.fd.fx-query.exact); 

%display_sol(obj, query);

figure(1);
title('Graph for arctanx, actual, forward difference');
hold on;
%axis([start-1 end+1 -inf inf]);
plot(query.v, query.exact, 'color', 'c', 'LineWidth', 5);
plot(query.v, query.fd.fx, 'color', 'm', 'LineWidth', 3);
hold off;

figure(2);
plot(query.v, query.diff);
title('Absolute differnce graph');
disp('Coefficients');


sym x;
input = 1:0.5:6;
output = atan(input);
n = size(input, 2);

fd_table = forward_difference_table2( output, n);
fx = 0;
syms x;
mu = (x-input(1))/(input(2)-input(1)); 
for i =1:n
    fx= fx+fd_table(1, i)*ncr_for2(mu, i-1); 
end

coeff = sym2poly(fx);
coeff = coeff';
disp(coeff);

function fd = newtons_forward_difference(obj, input)
    s = size(input, 2);
    fd.fd_table = forward_difference_table(obj);
    fd.fx = zeros(1, s);
    mu = (input-obj.x(1))/(obj.x(2)-obj.x(1));
    for i =1:obj.n
        fd.fx= fd.fx+fd.fd_table(1, i)*ncr_for(mu, i-1);
    end
end



function fd_table = forward_difference_table(obj)
    fd_table = zeros(obj.n, obj.n);
    fd_table(:, 1) = obj.y';
    for j = 2 : obj.n
        for i = 1 : (obj.n - j+1)
            fd_table(i,j) = fd_table(i + 1, j -1) - fd_table(i, j - 1);
        end
    end
end


function f = ncr_for(u, r)
    s = size(u, 2);
    f=ones(1, s);
    for i=0:r-1
        f = (f.*(u-i))/(i+1);
    end
end


function display_sol(obj, query)
    disp(obj);
    disp('Calculating for:');
    disp(query.v);
    disp('Exact Answer');
    disp(query.exact);
    disp('Solution');
    disp(query.fd.fx);
    
end

%coef

function fd_table = forward_difference_table2(output, n)
    fd_table = zeros(n, n);
    fd_table(:, 1) = output';
    for j = 2 : n
        for i = 1 : (n - j+1)
            fd_table(i,j) = fd_table(i + 1, j -1) - fd_table(i, j - 1);
        end
    end
end


function f = ncr_for2(u, r)
    syms x;
    f=1;
    for i=0:r-1
        f = (f*(u-i))/(i+1);
    end
end






