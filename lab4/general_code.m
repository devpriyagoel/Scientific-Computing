clear all
clc
clf
close all
disp('define nodal points for exp(x)');
start = input('Enter start value ');
increment = input('Enter increment ');
endi = input('Enter end value ');


obj.x = start:increment:endi;
obj.f = @(x) exp(x);
obj.y = obj.f(obj.x);
obj.n = size(obj.x, 2);

disp('define query array');
start = input('Enter start value ');
increment = input('Enter increment ');
endi = input('Enter end value ');

query.v = start:increment:endi;
query.fd = newtons_forward_difference(obj, query.v);
query.bd = newtons_backward_difference(obj, query.v);
query.exact = obj.f(query.v);

display_sol(obj, query);

figure(1);
title('Graph for exp(x), actual, forward difference');
hold on;
%axis([start-1 end+1 -inf inf]);
plot(query.v, query.exact, 'color', 'c', 'LineWidth', 5);
plot(query.v, query.fd.fx, 'color', 'm', 'LineWidth', 3);
%plot(query.v, query.bd.fx, 'color', 'k', 'Linewidth', 1);
hold off;

disp('define nodal points for 1/1+x^2');
start = input('Enter start value ');
increment = input('Enter increment ');
endi = input('Enter end value ');


obj2.x = start:increment:endi;
obj2.f = @(x) 1./(1+x.^2);
obj2.y = obj2.f(obj2.x);
obj2.n = size(obj2.x, 2);

disp('define query array');
start = input('Enter start value ');
increment = input('Enter increment ');
endi = input('Enter end value ');

query2.v = start:increment:endi;
query2.fd = newtons_forward_difference(obj2, query2.v);
query2.bd = newtons_backward_difference(obj2, query2.v);
query2.exact = obj2.f(query2.v);

display_sol(obj2, query2);

figure(2);
title('Graph for 1/1+x^2, actual, forward difference');
hold on;
%axis([start-1 end+1 -inf inf]);
plot(query2.v, query2.exact, 'color', 'c', 'LineWidth', 5);
plot(query2.v, query2.fd.fx, 'color', 'm', 'LineWidth', 3);
%plot(query.v, query.bd.fx, 'color', 'k', 'Linewidth', 1);
hold off;

input('press enter to exit');

function fd = newtons_forward_difference(obj, input)
    s = size(input, 2);
    fd.fd_table = forward_difference_table(obj);
    fd.fx = zeros(1, s);
    mu = (input-obj.x(1))/(obj.x(2)-obj.x(1));
    for i =1:obj.n
        fd.fx= fd.fx+fd.fd_table(1, i)*ncr_for(mu, i-1); 
    end
end

function bd = newtons_backward_difference(obj, input)
    s = size(input, 2);
    bd.bd_table = backward_difference_table(obj);
    bd.fx = zeros(1, s);
    mu = (input-obj.x(obj.n))/(obj.x(2)-obj.x(1));
    for i =1:obj.n
        bd.fx= bd.fx+bd.bd_table(obj.n, i)*ncr_back(mu, i-1); 
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

function bd_table = backward_difference_table(obj)
    bd_table = zeros(obj.n, obj.n);
    bd_table(:, 1) = obj.y';
    for j = 2:obj.n
        for i = obj.n:-1:j 
            bd_table(i,j) = bd_table(i, j -1) - bd_table(i-1, j-1);
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

function f = ncr_back(u, r)
    s = size(u, 2);
    f=ones(1, s);
    for i=0:r-1
        f = (f.*(u+i))/(i+1);
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
    disp('Solution');
    disp(query.bd.fx);
    
end
