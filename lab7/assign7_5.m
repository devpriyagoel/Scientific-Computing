clc;
clear all;
close all;

a = @(x) (x.^2 .* log(x));
syms x;
aa = (x.^2 .* log(x));

fprintf('\nFor a.\n');
fprintf('\nDegree\tGaussian quadrature\t\tExact\t\tError\n\n');
A = cell(1, 4);
for i=2:5
    A{i-1} = Gaussian_quadrature(a,aa,1,1.5,i);     
    fprintf('%d\t\t%0.10f\t\t%0.10f\t%.15f\n',i,A{i-1}.approx,A{i-1}.actual,A{i-1}.error);
end

b = @(x) ((exp(3.*x)).*(sin(2.*x)));
syms x;
bb = ((exp(3.*x)).*(sin(2.*x)));

fprintf('\nFor b.\n');
fprintf('\nDegree\tGaussian quadrature\t\tExact\t\tError\n\n');
B = cell(1, 4);
for i=2:5
    B{i-1} = Gaussian_quadrature(b,bb,0,pi/4,i);
    fprintf('%d\t\t%0.10f\t\t%0.10f\t%.15f\n',i,B{i-1}.approx,B{i-1}.actual,B{i-1}.error);
end

c = @(x) (2./(x.^2-4));
syms x;
cc = (2./(x.^2-4));

fprintf('\nFor c.\n');
fprintf('\nDegree\tGaussian quadrature\t\tExact\t\tError\n\n');

C = cell(1, 4);
for i=2:5
    C{i-1} = Gaussian_quadrature(c,cc,0,0.35,i);
    fprintf('%d\t\t%0.10f\t\t%0.10f\t%.15f\n',i,C{i-1}.approx,C{i-1}.actual,C{i-1}.error);
end

d = @(x) ((2.*x)./(x.^2-4));
syms x;
dd = ((2.*x)./(x.^2-4));

fprintf('\nFor d.\n');
fprintf('\nDegree\tGaussian quadrature\t\tExact\t\tError\n\n')

D = cell(1, 4);

for i=2:5
    D{i-1} = Gaussian_quadrature(d,dd,1,1.6,i);
    fprintf('%d\t\t%0.10f\t\t%0.10f\t%.15f\n',i,D{i-1}.approx,D{i-1}.actual,D{i-1}.error);
end

function sol = Gaussian_quadrature(f,ff,a,b,degree)
    syms x
    roots =  vpasolve(legendreP(degree,x) == 0);
    roots = double(roots);    
    for l = 1:degree
        temp  = 1;
        for m = 1:degree
            if m~=l
                temp = temp*(x - roots(m))/(roots(l) - roots(m));
            end
        end
        c(l) = double(int(temp,-1,1));
        c = c';
    end 
    sol.approx = 0;
    for i = 1:degree
        x = 0.5*(((b-a)*roots(i))+(b+a));
        sol.approx = sol.approx + c(i).*f(x);
    end
    sol.approx = (b-a)./2.*(sol.approx);
    sol.actual = double(int(ff,a,b));
    sol.error = abs(sol.approx-sol.actual);
end


