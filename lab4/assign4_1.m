clear all
clc
clf
close all

obj.x = [1.0 1.5 2.0 2.5];
obj.y = [2.7183 4.4817 7.3819 12.1825];
obj.n = size(obj.x, 2);
obj.f = @(x)exp(x);

query.v = 2.25;
query.fd = newtons_forward_difference(obj, query.v);
query.bd = newtons_backward_difference(obj, query.v);
query.exact = obj.f(query.v);

display_sol(obj, query);

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
    disp('Forward Difference Table');
    disp(query.fd.fd_table);
    disp('Solution');
    disp(query.fd.fx);
    disp('Backward Difference Table');
    disp(query.bd.bd_table);
    disp('Solution');
    disp(query.bd.fx);
    
end
