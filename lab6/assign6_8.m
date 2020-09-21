clc;
clear all;
close all;

partition = [1 5/4 3/2 7/4 2];
f = [10 8 7 6 5];
trapezoid(partition,f);

function trapezoid(a,f)
    n = length(a);
    y = 1/2.*(a(2:n)-a(1:n-1)).*(f(1:n-1)+f(2:n));
    fprintf('The approximate value of the integral using using the composite trapezoidal rule is : %0.3f\n',sum(y));
end
