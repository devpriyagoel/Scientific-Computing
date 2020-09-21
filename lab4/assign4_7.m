clear all
clc
clf
close all

obj.x = [2 3 5];
obj.y = [1.5713 1.5719 1.5738];
obj.n = size(obj.x, 2);

query.v = 4;
query.lm = lagrange_interpolation_method(obj, query.v);

disp('using second degree polynomial');
display_sol(obj, query);


obj2.x = [obj.x 6];
obj2.y = [obj.y 1.5751];
obj2.n = size(obj2.x, 2);

query2.v = 4;
query2.lm = lagrange_interpolation_method(obj2, query2.v);

disp('using third degree polynomial');
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

function display_sol(obj, query)
    disp(obj);
    disp('Calculating for:');
    disp(query.v);
    disp('Lagrange Interpolation');
    disp(query.lm);
end