clear all
clc
clf
close all

obj.x = 0.1:0.1:0.4;
obj.y = [-0.29004986 -0.56079734 -0.81401972 -1.0526302];
obj.n = size(obj.x, 2);

query.v = 0.18;
query.lm = lagrange_interpolation_method(obj, query.v);

disp('part a: ');
display_sol(obj, query);


obj2.x = -1:0.5:0.5;
obj2.y = [0.86199480 0.95802009 1.0986123 1.2943767];
obj2.n = size(obj2.x, 2);

query2.v = 0.25;
query2.lm = lagrange_interpolation_method(obj2, query2.v);

disp('part b');
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