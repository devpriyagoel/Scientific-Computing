clc;
clear all;
close all;

%problem A

f = @(x, y)0;
eqx1 = @(x)0;
eqx2 = @(x)x;
eqy1 = @(y)0;
eqy2 = @(y)y;
a=0;
b=1;
h = 2e-1;
ex = @(x, y)x*y;

global fig;
fig = 1;

error_tableg = zeros(5, 1);
error_tablej = zeros(5, 1);
n = zeros(5, 1);
for i= 1:5
    xd = a:h:b;
    yd = a:h:b;
    nx = size(xd, 2);
    ny = size(yd, 2);
    n(i) = h;
    [error_tableg(i), error_tablej(i)] = five_point_stencil_gaussjacobi(h, xd, yd, nx, ny, f, ex, eqx1, eqx2, eqy1, eqy2);
    h = h./2;
end
figure;
loglog(n,error_tableg,'LineWidth',1.5);
xlabel('\Delta x(=\Delta y)');
ylabel('Max Absolute Error');
title('loglog plot: \Delta x(=\Delta y) vs max error, 5-point stencil (Gauss Siedel) for problem (a)');
saveas(gcf,sprintf('q2_af%d.png',fig));
fig = fig+1;

figure;
loglog(n,error_tablej,'LineWidth',1.5);
xlabel('\Delta x(=\Delta y)');
ylabel('Max Absolute Error');
title('loglog plot: \Delta x(=\Delta y) vs max error, 5-point stencil (Jacobi) for problem (a)');
saveas(gcf,sprintf('q2_af%d.png',fig));
fig = fig+1;
h=2e-2;
for i= 1:1
    xd = a:h:b;
    yd = a:h:b;
    nx = size(xd, 2);
    ny = size(yd, 2);
    str = strcat('problem(a) h = k = ', num2str(h));
    five_point_stencil_gaussjacobif(h, xd, yd, nx, ny, f, ex, eqx1, eqx2, eqy1, eqy2, str);
end

function [max_errorg, max_errorj] = five_point_stencil_gaussjacobi(h, xd, yd, nx, ny, f, ex, eqx1, eqx2, eqy1, eqy2)
    B = zeros(nx*ny, 1);
    tabg = zeros(nx, ny);
    tabj1 = zeros(nx, ny);
    tabj2 = zeros(nx, ny);
    tabex = zeros(nx, ny);
    errorg = zeros(nx, ny);
    errorj = zeros(nx, ny);
    for i=1:ny  
        tabg(1, i) = eqy1(yd(i));  
        tabj1(1, i) = eqy1(yd(i));  
        tabj2(1, i) = eqy1(yd(i)); 
        B(i) = eqy1(yd(i)); 
        tabg(nx, i) = eqy2(yd(i));
        tabj1(nx, i) = eqy2(yd(i));
        tabj2(1, i) = eqy1(yd(i)); 
        B((nx-1)*ny+i) = eqy2(yd(i));
    end

    for i=1:nx 
        tabg(i, 1) = eqx1(xd(i));
        tabj1(i, 1) = eqx1(xd(i));
        tabj2(1, i) = eqy1(yd(i)); 
        B((i-1)*ny+1) = eqx1(xd(i));
        tabg(i, ny) = eqx2(xd(i));
        tabj1(i, ny) = eqx2(xd(i));
        tabj2(1, i) = eqy1(yd(i)); 
        B((i)*ny) = eqx2(xd(i));
    end

    for i=2:nx-1
        for j=2:ny-1
            B((i-1)*ny+j) = h*h*f(xd(i), yd(j));
        end
    end
    for i=1:nx
        for j=1:ny
            tabex(i, j) = ex(xd(i), yd(j));
        end
    end
    it = 0;
    while(it<1e5)
        for i=2:nx-1
            for j=2:ny-1
                tabg(i, j) = (-0.25)*(B((i-1)*ny+j)-tabg(i-1, j)-tabg(i, j-1)-tabg(i, j+1)-tabg(i+1, j));
                tabj2(i, j) = (-0.25)*(B((i-1)*ny+j)-tabj1(i-1, j)-tabj1(i, j-1)-tabj1(i, j+1)-tabj1(i+1, j));
            end
        end
        for i=2:nx-1
            for j=2:ny-1
                tabj1(i, j) = tabj2(i, j);
            end
        end
        it =it+1;
    end
    max_errorg=0;
    max_errorj=0;
    for i=1:nx
        for j=1:ny
            errorg(i, j) = abs(tabg(i, j)-tabex(i, j));
            errorj(i, j) = abs(tabj1(i, j)-tabex(i, j));
            max_errorg = max([max_errorg errorg(i, j)]);
            max_errorj = max([max_errorj errorj(i, j)]);
        end
    end
end


function [max_errorg, max_errorj] = five_point_stencil_gaussjacobif(h, xd, yd, nx, ny, f, ex, eqx1, eqx2, eqy1, eqy2, str)
    global fig;
    B = zeros(nx*ny, 1);
    tabg = zeros(nx, ny);
    tabj1 = zeros(nx, ny);
    tabj2 = zeros(nx, ny);
    tabex = zeros(nx, ny);
    errorg = zeros(nx, ny);
    errorj = zeros(nx, ny);
    for i=1:ny  
        tabg(1, i) = eqy1(yd(i));  
        tabj1(1, i) = eqy1(yd(i));  
        tabj2(1, i) = eqy1(yd(i)); 
        B(i) = eqy1(yd(i)); 
        tabg(nx, i) = eqy2(yd(i));
        tabj1(nx, i) = eqy2(yd(i));
        tabj2(1, i) = eqy1(yd(i)); 
        B((nx-1)*ny+i) = eqy2(yd(i));
    end

    for i=1:nx 
        tabg(i, 1) = eqx1(xd(i));
        tabj1(i, 1) = eqx1(xd(i));
        tabj2(1, i) = eqy1(yd(i)); 
        B((i-1)*ny+1) = eqx1(xd(i));
        tabg(i, ny) = eqx2(xd(i));
        tabj1(i, ny) = eqx2(xd(i));
        tabj2(1, i) = eqy1(yd(i)); 
        B((i)*ny) = eqx2(xd(i));
    end

    for i=2:nx-1
        for j=2:ny-1
            B((i-1)*ny+j) = h*h*f(xd(i), yd(j));
        end
    end
    for i=1:nx
        for j=1:ny
            tabex(i, j) = ex(xd(i), yd(j));
        end
    end
    it = 0;
    while(it<1e5)
        for i=2:nx-1
            for j=2:ny-1
                tabg(i, j) = (-0.25)*(B((i-1)*ny+j)-tabg(i-1, j)-tabg(i, j-1)-tabg(i, j+1)-tabg(i+1, j));
                tabj2(i, j) = (-0.25)*(B((i-1)*ny+j)-tabj1(i-1, j)-tabj1(i, j-1)-tabj1(i, j+1)-tabj1(i+1, j));
            end
        end
        for i=2:nx-1
            for j=2:ny-1
                tabj1(i, j) = tabj2(i, j);
            end
        end
        it =it+1;
    end
    max_errorg=0;
    max_errorj=0;
    for i=1:nx
        for j=1:ny
            errorg(i, j) = abs(tabg(i, j)-tabex(i, j));
            errorj(i, j) = abs(tabj1(i, j)-tabex(i, j));
            max_errorg = max([max_errorg errorg(i, j)]);
            max_errorj = max([max_errorj errorj(i, j)]);
        end
    end
    [X,Y]=meshgrid(xd,yd);
    figure;
    surf(X,Y,tabex);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Surface Plot: Exact solution for ', str));
    saveas(gcf,sprintf('q2_af%d.png',fig));
    fig = fig + 1;
    figure;
    surf(X,Y,tabg);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Solution by 5-point stencil(Gauss Seidel) for ', str));
    saveas(gcf,sprintf('q2_af%d.png',fig));
    fig = fig + 1;
    figure;
    surf(X,Y,tabj1);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Solution by 5-point stencil (Jacobi) for ', str));
    saveas(gcf,sprintf('q2_af%d.png',fig));
    fig = fig + 1;
    figure;
    contour(X,Y,tabex);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Contour Plot: Exact solution for ', str));
    saveas(gcf,sprintf('q2_af%d.png',fig));
    fig = fig + 1;
    figure;
    contour(X,Y,tabg);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Solution by 5-point stencil(Gauss Seidel) for ', str));
    saveas(gcf,sprintf('q2_af%d.png',fig));
    fig = fig + 1;
    figure;
    contour(X,Y,tabj1);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Solution by 5-point stencil (Jacobi) for ', str));
    saveas(gcf,sprintf('q2_af%d.png',fig));
    fig = fig + 1;
    figure;
    surf(X,Y,errorg);
    xlabel('x');
    ylabel('y');
    zlabel('Absolute Error');
    title(strcat('Error: 5-point stencil(Gauss Seidel) for ', str));
    saveas(gcf,sprintf('q2_af%d.png',fig));
    fig = fig + 1;
    figure;
    surf(X,Y,errorj);
    xlabel('x');
    ylabel('y');
    zlabel('Absolute Error');
    title(strcat('Error: 5-point stencil (Jacobi) for ', str));
    saveas(gcf,sprintf('q2_af%d.png',fig));
    fig = fig + 1;
end
