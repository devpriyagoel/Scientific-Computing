clc;
clear all;
close all;

Q = cell(4, 1);
Q{1}.f = @(x)(x.^2.*exp(-1.*x.*x));
Q{2}.f = @(x)(1./x./log(x));
Q{3}.f = @(x)((x.^2).*log(x.^2+1));
Q{4}.f = @(x)((sin(x)).^2 - 2.*x.*sin(x) + 1);
Q{1}.a = 1;
Q{1}.b = 1.6;
Q{2}.a = 0;
Q{2}.b = pi/4;
Q{3}.a = .75;
Q{3}.b = 1.3;
Q{4}.a = exp(1);
Q{4}.b = exp(1)+1;

Q{1}.n = 8;
Q{1}.a = 0;
Q{1}.b = 2;
Q{1}.h = (Q{1}.b-Q{1}.a)/Q{1}.n;


Q{2}.n = 8;
Q{2}.a = exp(1);
Q{2}.b = exp(1)+2;
Q{2}.h = (Q{2}.b-Q{2}.a)/Q{2}.n;


Q{3}.n = 8;
Q{3}.a = 0;
Q{3}.b = 2;
Q{3}.h = (Q{3}.b-Q{3}.a)/Q{3}.n;


Q{4}.n = 8;
Q{4}.a = 0.75;
Q{4}.b = 1.75;
Q{4}.h = (Q{4}.b-Q{4}.a)/Q{4}.n;


for i=1:4
    Q{i}.trap = trapezoid(Q{i}.a, Q{i}.b, Q{i}.h, Q{i}.f);
    Q{i}.midp = midpoint(Q{i}.a, Q{i}.b, Q{i}.h, Q{i}.f);
    Q{i}.simp = simpson(Q{i}.a, Q{i}.b, Q{i}.h, Q{i}.f);
end
for i=1:4
    fprintf('Value of integral using composite trapezoidal: %f\n',Q{i}.trap);
    fprintf('Value of integral using composite midpoint: %f\n',Q{i}.midp);
    fprintf('Value of integral using composite simpsons: %f\n',Q{i}.simp);
end
function solution = trapezoid(a,b,h,f)
    partition = a:h:b;
    f_at_x = sum(f(partition)) + sum(f(partition(2:(length(partition)-1))));
    solution = h/2*f_at_x;
end

function solution = midpoint(a,b,h,f)
    partition = a+h/2:h:b-h/2;
    f_at_x = sum(f(partition));
    solution = h*f_at_x;
end

function solution = simpson(a,b,h,f)
    partition = a+h/2:h:b-h/2;
    f_at_x = sum(f(partition-h/2) + 4.*f(partition) + f(partition+h/2));
    solution = h/6*f_at_x;
end

