clear all
clc
clf
close all

obj.x = [1950 1960 1970 1980 1990 2000];
obj.y = [151326 179323 203302 226542 249633 281422];
obj.n = size(obj.x, 2);

query.v = [1940 1975 2020];
query.dd = divided_difference_method(obj, query.v);
display_sol(obj, query);


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
    disp('Divided Difference Table');
    disp(query.dd.dd_table);
    disp('Divided Difference Interpolation');
    disp(query.dd.fx);
end