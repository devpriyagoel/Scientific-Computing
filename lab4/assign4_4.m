clear all
clc
clf
close all

obj.x = [0.0 0.1 0.3 0.6 1.0];
obj.y = [-6.00000 -5.89483 -5.65014 -5.17788 -4.28172];
obj.n = size(obj.x, 2);

query.v = 0.2;
query.lm = lagrange_interpolation_method(obj, query.v);
query.dd = divided_difference_method(obj, query.v);

display_sol(obj, query);

% adding f(1.1) =−3.99583

obj2.x = [obj.x 1.1];
obj2.y = [obj.y -3.99583];
obj2.n = size(obj2.x, 2);

query2.v = 0.2;
query2.lm = lagrange_interpolation_method(obj2, query2.v);
query2.dd = divided_difference_method(obj2, query2.v);

display_sol(obj2, query2);

function fx = lagrange_interpolation_method(obj, input)
   s = size(input, 2);
   fx = zeros(1, s);
   for i=1:obj.n
       fx = fx+obj.y(i)*(li_term(obj, input, i));
   end
end


function li = li_term(obj, input, t)        
    s = size(input, 2);
    li=ones(1, s);
    for i=1:obj.n
        if(t==i)
            continue;
        end
        li = li.*(input-obj.x(i));
    end
    for i=1:obj.n
        if(t==i)
            continue;
        end
        li = li/(obj.x(t)-obj.x(i));
    end
end

function dd =  divided_difference_method(obj, input)
    dd.dd_table = divided_difference_table(obj);
    s = size(input, 2);
    dd.fx=zeros(1, s);
    prod = ones(1, s);
    for i=1:obj.n
        dd.fx = dd.fx + dd.dd_table(1, i).*prod;
        prod = prod.*(input-obj.x(i));
    end
end

function dd_table = divided_difference_table(obj)
    dd_table = zeros(obj.n, obj.n);
    dd_table(:, 1) = obj.y';
    for j = 2 : obj.n
        for i = 1 : (obj.n - j + 1)
            dd_table(i,j) = (dd_table(i + 1, j - 1) - dd_table(i, j - 1)) / (obj.x(i + j - 1) - obj.x(i));
        end
    end
end
function display_sol(obj, query)
    disp(obj);
    disp('Calculating for:');
    disp(query.v);
    disp('Lagrange Interpolation');
    disp(query.lm);
    disp('divided difference table');
    disp(query.dd.dd_table);
    disp('Divided Difference Interpolation');
    disp(query.dd.fx);
end