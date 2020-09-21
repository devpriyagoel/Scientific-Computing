clc;
clear all;
close all;
a = -4;
b = 4;
obj.f = @(x) (log(exp(x)+2));
obj.x = [-1,-0.5,0,0.5]; 
obj.y = [0.86199480,0.95802009,1.0986123,1.2943767];
obj.z = [0.15536240,0.23269654,0.33333333,0.45186776];
obj.dd_table = divided_difference(obj.x,obj.y,obj.z);
hermite = Hermite(obj.dd_table);

x = 0.25;
fx_hermite = hermite(x);
fx_actual = obj.f(x);
error = abs(fx_actual-fx_hermite);
fprintf('Using Hermite interpolation at x = %f : \nf(x) = %e\nActual value = %e\nAbsolute error = %e\n',x,fx_hermite,fx_actual,error);

hold on;
fplot(@(x)obj.f(x),[a,b],'color','magenta','LineWidth',1);
fplot(@(x)hermite(x),[a,b],'color','blue');
legend('Actual function','Hermite interpolation polynomial')
hold off;

function X = divided_difference(x,y,z)
    n = 2*length(x);
    for i = 1:n
        if(mod(i,2)==1)
            X(i,1)=x(ceil(i/2));
            X(i,2)=y(ceil(i/2));
            
        else
            if (i~=2)
                X(i-1,3) = (X(i-1,2)-X(i-2,2))/ (X(i-1,1)-X(i-2,1));
            end
            X(i,1)=x(i/2);
            X(i,2)=y(i/2);
            X(i,3)=z(i/2);
        end
    end
    for j = 4:n+1
        for i = j-1:n
            X(i,j) = (X(i,j-1) - X(i-1,j-1))/(X(i,1) - X(i-j+2,1));
        end
    end
end

function hermite = Hermite(X)
    n = size(X, 1);
    syms x
    for i = 2:n+1
        Y(i-1) = X(1,i);
    end
    hermite = Y(1);
    for i = 1:n-1
        temp = 1;
        for j = 1:i
            temp = temp*(x-X(j,1));
        end
        hermite = hermite+temp*X(i+1,i+2);
    end
    hermite = matlabFunction(hermite);
end

