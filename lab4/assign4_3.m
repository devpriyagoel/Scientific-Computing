clear all
clc
clf
close all

obj.x = 1:0.1:1.4;
obj.y = log(obj.x);
obj.n = size(obj.x, 2);

query.v = 1:0.01:1.4;
query.lm = lagrange_interpolation_method(obj, query.v);
query.diff = abs(log(query.v)-query.lm);

figure(1);

plot(query.v, query.diff);
title('Absolute difference graph');

display_sol(obj, query);
disp(max(query.diff));
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