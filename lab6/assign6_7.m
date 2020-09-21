clc;
clear all;
close all;

obj1.f = @(x)(4./(1+x.^2));
obj2.f = @(x) (((1-x.^2).^(0.5))-x);

obj1.a = 0;
obj1.b = 1;
obj1.ans = simpson(obj1.a,obj1.b,obj1.f,1);

obj2.a = 0;
obj2.b = 1/(2)^0.5;
obj2.ans = simpson(obj2.a,obj2.b,obj2.f,2);

printSol(obj1.ans);
printSol(obj2.ans);
function X = simpson(a,b,f,i)
    n = 1;
    actual = integral(f,a,b);
    approx = (b-a)/6*(f(a)+4*f((a+b)/2)+f(b));
    while (abs(actual-approx)>=0.5*1e-10)
        X(n,1)=n;X(n,2)=actual;X(n,3)=approx;X(n,4)=abs(actual-approx);
        n = n+1;
        h = (b-a)/n;
        part = [a+h/2:h:b-h/2];
        f_at_x = sum(f(part-h/2) + 4.*f(part) + f(part+h/2));
        approx = h/6*f_at_x;
    end
    figure(i);
    plot(X(:,1),log(X(:,4)));
    title('Number of iterations vs log(Error)');
    xlabel('n');
    ylabel('log(Error)');
    X(n,1)= n; X(n, 2) = actual;X(n,3)=approx;X(n,4)=actual-approx;
end

function printSol(X)
    fprintf('\nn\tActual\t\t\tSimpson for n partitions\tError\n');
    n = size(X);
    for i=1:n
        fprintf('\n%d\t%0.12f\t\t%0.13f\t\t\t%0.13f\n',X(i, 1),X(i, 2),X(i, 3),X(i, 4));
    end
    
end