clear all;
close all;
clc;
%heat conductivity
q1.f = @(ft, fx) ft-(4/(pi*pi))*fx;
q1.eqx = @(x)sin((pi*x)/4).*(1+2*cos((pi*x)/4));
q1.xa = 0;
q1.xb = 4;
q1.ta = 0;
q1.tb = 1;
q1.eqt1 = @(t)t-t;
q1.eqt2 = @(t)t-t;
q1.ex = @(x, t)exp(-t).*(sin((pi*x)/2))+exp(-t/4).*sin((pi*x)/4);
s1 = 4/(pi*pi);
f = solve(q1, 1, s1);


q2.f = @(ft, fx) ft-fx;
q2.eqx = @(x)sin(pi*x);
q2.xa = 0;
q2.xb = 1;
q2.ta = 0;
q2.tb = 1;
q2.eqt1 = @(t)t-t;
q2.eqt2 = @(t)t-t;
q2.ex = @(x, t)exp(-pi*pi*t).*(sin(pi*x));
s2=1;
f = solve(q2, f, s2);

q3.f = @(ft, fx) ft-fx;
q3.eqx = @(x)sin((pi*x)/2)+0.5*sin(2*pi*x);
q3.xa = 0;
q3.xb = 1;
q3.ta = 0;
q3.tb = 1;
q3.eqt1 = @(t)t-t;
q3.eqt2 = @(t)exp(-pi*pi*t/4);
q3.ex = @(x, t)exp(-pi*pi*t*0.25).*sin(pi*x*0.5)+0.5*exp(-4*pi*pi*t)*sin(2*pi*x);
s3=1;
solve(q3, f, s3);

function f = solve(q, f, s)
    h=1;
    errorftcs = zeros(1, 4);
    errorbtcs = zeros(1, 4);
    errorcrank = zeros(1, 4);
    errorrich = zeros(1, 4);
    errorfrank = zeros(1, 4);
    n = zeros(1, 4);
    w=1;
    while(h>0.1)
        k=1e-2;
        td = q.ta:k:q.tb;
        xd = q.xa:h:q.xb;
        nt = size(td, 2);
        nx = size(xd, 2);
        u = zeros(nx, nt);
        ex= zeros(nx, nt);
        error= zeros(nx, nt);
        lambda = (k/(h*h))*s;
        for i=1:nx
            u(i, 1) = q.eqx(xd(i));
        end
        for i=1:nt
            u(1, i) = q.eqt1(td(i));
            u(nx, i) = q.eqt2(td(i));
        end
        for i=1:nt
            for j=1:nx
                ex(j, i) = q.ex(xd(j), td(i));
            end
        end
        errorftcs(w) = ftcs(nx, nt, xd, td, lambda, u, ex, error, f);
        errorbtcs(w) = btcs(nx, nt, xd, td, lambda, u, ex, error, f+1);
        errorcrank(w) = crank(nx, nt, xd, td, lambda, u, ex, error, f+2);
        errorrich(w) = richardson(nx, nt, xd, td, lambda, u, ex, error, f+3);
        errorfrank(w) = frankel(nx, nt, xd, td, lambda, u, ex, error, f+4);    
        n(w) = nx;
        f=f+5;
        h = h/2;
        w=w+1;
    end
    figure(f);
    hold on;
    plot(log(n), log(errorftcs));
    plot(log(n), log(errorbtcs));
    plot(log(n), log(errorcrank));
    plot(log(n), log(errorrich));
    plot(log(n), log(errorfrank));
    title('loglog error graph for all schemes');
    legend('FTCS', 'BTCS', 'Crank-Nicolson', 'Richardson', 'DuFort-Frankel');
    hold off;
    f=f+1;
end

function max_error = ftcs(nx, nt, xd, td, lambda, u, ex, error, f)
    for i=2:nt
        for j=2:nx-1
            u(j, i) = lambda*u(j-1, i-1)+(1-2*lambda)*u(j, i-1)+lambda*u(j+1, i-1);
        end
    end
    for i=1:nt
        for j=1:nx
            error(j, i) = abs(ex(j, i)-u(j, i));
        end
    end
    max_error = max(error(:, nt));
    [X, T] = meshgrid(td, xd);
    figure(f);
    subplot(2, 2, 1);
    surf(X, T, u);
    title(['FTCS: Approximate solution lambda = ', num2str(lambda)]);
    subplot(2, 2, 2);
    surf(X, T, ex);
    title('Exact solution');
    subplot(2, 2, 3);
    surf(X, T, error);
    title('Absolute Error');
    subplot(2, 2, 4);
    hold on;
    title('Graph at time t Approx vs Exact');
    plot(xd(:), u(:, nt));
    plot(xd(:), ex(:, nt));
    legend('Approximate', 'Exact');
    hold off;
end
    
function max_error = btcs(nx, nt, xd, td, lambda, u, ex, error, f)
    A = zeros(nx, nx);
    A(1, 1)=1;
    A(nx, nx)=1;
    for j=2:nx-1
        A(j, j) = 1+lambda;
        A(j-1, j) = -lambda/2;
        A(j+1, j) = -lambda/2;
    end
    for i=2:nt
        B = zeros(nx);
        B(1) = u(1, i);
        B(nx) =  u(nx, i);
        for j=2:nx-1
            B(j) = u(j, i-1);
        end
        
        C = A\B;
        for j=2:nx-1
            u(j, i) = C(j);
        end
    end
    for i=1:nt
        for j=1:nx
            error(j, i) = abs(ex(j, i)-u(j, i));
        end
    end
    max_error = max(error(:, nt));
    [X, T] = meshgrid(td, xd);
    figure(f);
    subplot(2, 2, 1);
    surf(X, T, u);
    title(['BTCS: Approximate solution lambda = ', num2str(lambda)]);
    subplot(2, 2, 2);
    surf(X, T, ex);
    title('Exact solution');
    subplot(2, 2, 3);
    surf(X, T, error);
    title('Absolute Error');
    subplot(2, 2, 4);
    hold on;
    title('Graph at time t Approx vs Exact');
    plot(xd(:), u(:, nt));
    plot(xd(:), ex(:, nt));
    legend('Approximate', 'Exact');
    hold off;
end

function max_error = crank(nx, nt, xd, td, lambda, u, ex, error, f)
    A = zeros(nx, nx);
    A(1, 1)=1;
    A(nx, nx)=1;
    for j=2:nx-1
        A(j, j) = 1+lambda;
        A(j-1, j) = -lambda/2;
        A(j+1, j) = -lambda/2;
    end
    for i=2:nt

        B = zeros(nx);
        B(1) = u(1, i);
        B(nx) =  u(nx, i);
        for j=2:nx-1
            B(j) = u(j, i-1);
        end
        
        C = A\B;
        for j=2:nx-1
            u(j, i) = C(j);
        end
    end
    for i=1:nt
        for j=1:nx
            error(j, i) = abs(ex(j, i)-u(j, i));
        end
    end
    max_error = max(error(:, nt));
    [X, T] = meshgrid(td, xd);
    figure(f);
    subplot(2, 2, 1);
    surf(X, T, u);
    title(['Crank-Nicolson Approximate solution lambda = ', num2str(lambda)]);
    subplot(2, 2, 2);
    surf(X, T, ex);
    title('Exact solution');
    subplot(2, 2, 3);
    surf(X, T, error);
    title('Absolute Error');
    subplot(2, 2, 4);
    hold on;
    title('Graph at time t Approx vs Exact');
    plot(xd(:), u(:, nt));
    plot(xd(:), ex(:, nt));
    legend('Approximate', 'Exact');
    hold off;
end

function max_error = richardson(nx, nt, xd, td, lambda, u, ex, error, f)
    for j=2:nx-1
        u(j, 2) = lambda*u(j-1, 1)+(1-2*lambda)*u(j, 1)+lambda*u(j+1, 1);
    end
    for i=3:nt
        for j=2:nx-1
            u(j, i) = 2*lambda*(u(j+1, i-1)-2*u(j, i-1)+u(j-1, i-1))+u(j, i-2);
        end
    end
    for i=1:nt
        for j=1:nx
            error(j, i) = abs(ex(j, i)-u(j, i));
        end
    end
    max_error = max(error(:, nt));
    [X, T] = meshgrid(td, xd);
    figure(f);
    subplot(2, 2, 1);
    surf(X, T, u);
    title(['Richardson: Approximate solution lambda = ', num2str(lambda)]);
    subplot(2, 2, 2);
    surf(X, T, ex);
    title('Exact solution');
    subplot(2, 2, 3);
    surf(X, T, error);
    title('Absolute Error');
    subplot(2, 2, 4);
    hold on;
    title('Graph at time t Approx vs Exact');
    plot(xd(:), u(:, nt));
    plot(xd(:), ex(:, nt));
    legend('Approximate', 'Exact');
    hold off;
end
   
function max_error = frankel(nx, nt, xd, td, lambda, u, ex, error, f)
    for j=2:nx-1
        u(j, 2) = lambda*u(j-1, 1)+(1-2*lambda)*u(j, 1)+lambda*u(j+1, 1);
    end
    for i=3:nt
        for j=2:nx-1
            u(j, i) = ((1-2*lambda)/(1+2*lambda))*u(j, i-2)+((2*lambda)/(1+2*lambda))*(u(j+1, i-1)+u(j-1, i-1));
        end
    end
    for i=1:nt
        for j=1:nx
            error(j, i) = abs(ex(j, i)-u(j, i));
        end
    end
    max_error = max(error(:, nt));
    [X, T] = meshgrid(td, xd);
    figure(f);
    subplot(2, 2, 1);
    surf(X, T, u);
    title(['DuFort-Frankel: Approximate solution lambda = ', num2str(lambda)]);
    subplot(2, 2, 2);
    surf(X, T, ex);
    title('Exact solution');
    subplot(2, 2, 3);
    surf(X, T, error);
    title('Absolute Error');
    subplot(2, 2, 4);
    hold on;
    title('Graph at time t Approx vs Exact');
    plot(xd(:), u(:, nt));
    plot(xd(:), ex(:, nt));
    legend('Approximate', 'Exact');
    hold off;
end
   

