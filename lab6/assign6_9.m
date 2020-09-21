clc;
clear all;
close all;

f = @(x) (x.*log(x));

hold on;
Trapezoidal = trapezoidal(1,2,f);
Simpson = simpson(1,2,f);
Midpoint = midpoint(1,2,f);
legend('Trapezoidal','Simpson','Midpoint');
hold off;
disp('Numerical Approximation of the given Integral');
fprintf('\n\nFor Simpson Rule : \n\n');
printSol(Simpson);
fprintf('\n\nFor Trapezoidal Rule : \n\n');
printSol(Trapezoidal);
fprintf('\n\nFor Midpoint Rule : \n\n');
printSol(Midpoint);
function X = simpson(a,b,f)
    n = 1;
    actual = integral(f,a,b);
    approx = (b-a)/6*(f(a)+4*f((a+b)/2)+f(b));
    while (abs(actual-approx)>=1e-10)
        X(n,1)=n;X(n,2)=actual;X(n,3)=approx;X(n,4)=abs(actual-approx);
        n = n+1;
        h = (b-a)/n;
        part = [a+h/2:h:b-h/2];
        f_at_x = sum(f(part-h/2) + 4.*f(part) + f(part+h/2));
        approx = h/6*f_at_x;
    end
    figure(1);
    plot(X(:,1),log(X(:,4)),'color','m');
    title('log(Error) vs Number of iterations');
    xlabel('n');
    ylabel('log(Error)');
    X(n,1)=n; X(n, 2)=actual;X(n,3)=approx;X(n,4)=actual-approx;
end

function X = trapezoidal(a,b,f)
    n = 1;
    actual = integral(f,a,b);
    approx = ((b-a)/2)*(f(a)+f(b));
    while (abs(actual-approx)>=1e-5)
        X(n,1)=n;X(n,2)=actual;X(n,3)=approx;X(n,4)=abs(actual-approx);
        n = n+1;
        h = (b-a)/n;
        part = [a:h:b];
        part1 = part(1:length(part)-1);
        part2 = part(2:length(part));
        part = part2 - part1;
        f1 = f(part1);
        f2 = f(part2);
        f_ = f1 + f2;
        areas = 1/2.*part.*f_;
        approx = sum(areas);
    end
    figure(1);
    plot(X(:,1),log(X(:,4)),'color','b');
    title('log(Error) vs Number of iterations');
    xlabel('n');
    ylabel('log(Error)');
    X(n,1)=n; X(n, 2)=actual;X(n,3)=approx;X(n,4)=actual-approx;
end

function X = midpoint(a,b,f)
    n = 1;
    actual = integral(f,a,b);
    approx = (b-a)*(f((a+b)/2));
    while (abs(actual-approx)>=1e-5)
        X(n,1)=n;X(n,2)=actual;X(n,3)=approx;X(n,4)=abs(actual-approx);
        n = n+1;
        h = (b-a)/n;
        part = [(a+(h/2)):h:(b-(h/2))];
        areas = h.*(f(part));
        approx = sum(areas);
    end
    figure(1);
    plot(X(:,1),log(X(:,4)),'color','r');
    title('log(Error) vs Number of iterations');
    xlabel('n');
    ylabel('log(Error)');
    X(n,1)=n; X(n, 2)=actual;X(n,3)=approx;X(n,4)=actual-approx;
end
function printSol(X)
    fprintf('\nn\tActual\t\t\tApproximation for n partitions\tError\n');
    n = size(X);
    for i=1:n
        fprintf('\n%d\t%0.12f\t\t%0.13f\t\t\t%0.13f\n',X(i, 1),X(i, 2),X(i, 3),X(i, 4));
    end
    
end