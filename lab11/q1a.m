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

error_table = zeros(5, 1);
n = zeros(5, 1);
for i= 1:5
    xd = a:h:b;
    yd = a:h:b;
    nx = size(xd, 2);
    ny = size(yd, 2);
    n(i) = h;
    error_table(i) = five_point_stencil_matrix(h, xd, yd, nx, ny, f, ex, eqx1, eqx2, eqy1, eqy2);
    h = h./2;
end
figure;
loglog(n,error_table,'LineWidth',1.5);
xlabel('\Delta x(=\Delta y)');
ylabel('Max Absolute Error');
title('loglog plot: \Delta x(=\Delta y) vs max error by five-point stencil for problem (a)');
saveas(gcf,sprintf('q1_af%d.png',fig));
fig = fig+1;
h=2e-2;
for i= 1:1
    xd = a:h:b;
    yd = a:h:b;
    nx = size(xd, 2);
    ny = size(yd, 2);
    str = strcat('problem(a) h = k = ', num2str(h));
    five_point_stencil_matrixf(h, xd, yd, nx, ny, f, ex, eqx1, eqx2, eqy1, eqy2, str);
end
function max_error = five_point_stencil_matrixf(h, xd, yd, nx, ny, f, ex, eqx1, eqx2, eqy1, eqy2, str)
    global fig;
    A = zeros(nx*ny, nx*ny);
    B = zeros(nx*ny, 1);
    for i=1:ny-1   %at x=0%
        A(i, i)=1;  %u11, u12, u13%
        B(i) = eqy1(yd(i));  
    end

    for i=1:ny  %at x=1%
        A((nx-1)*ny+i,(nx-1)*ny+i) = 1; %u41 u42%
        B((nx-1)*ny+i) = eqy2(yd(i));
    end

    for i=2:nx  %at y=0%
        A((i-1)*ny+1, (i-1)*ny+1) = 1;%u11 u21 u31%
        B((i-1)*ny+1) = eqx1(xd(i));
    end

    for i=1:nx %at y=4%
        A((i)*ny, (i)*ny) = 1;%u14 u24 %
        B((i)*ny) = eqx2(xd(i));
    end

    for i=2:nx-1
        for j=2:ny-1
            A((i-1)*ny+j, (i-1)*ny+j) = -4;
            A((i-1)*ny+j, (i-1)*ny+j-1) = 1;
            A((i-1)*ny+j, (i-1)*ny+j+1) = 1;
            A((i-1)*ny+j, (i)*ny+j) = 1;
            A((i-1)*ny+j, (i-2)*ny+j) = 1;
            B((i-1)*ny+j) = h*h*f(xd(i), yd(j));
        end
    end

    u = A\B;
    tab = zeros(nx, ny);
    tabex = zeros(nx, ny);
    error = zeros(nx, ny);
    max_error=0;
    for i=1:nx
        for j=1:ny
            tab(i, j) = u((i-1)*ny+j);
            tabex(i, j) = ex(xd(i), yd(j));
            error(i, j) = abs(tab(i, j)-tabex(i, j));
            max_error = max([max_error error(i, j)]);
        end
    end
    [X,Y]=meshgrid(xd,yd);
    figure;
    surf(X,Y,tabex);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Surface Plot: Exact solution for ', str));
    saveas(gcf,sprintf('q1_af%d.png',fig));
    fig = fig + 1;
    figure;
    surf(X,Y,tab);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Surface Plot: Solution by five-point stencil for ', str));
    saveas(gcf,sprintf('q1_af%d.png',fig));
    fig = fig + 1;
    figure;
    contour(X,Y,tabex);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Contour Plot: Exact solution for ', str));
    saveas(gcf,sprintf('q1_af%d.png',fig));
    fig = fig + 1;
    figure;
    contour(X,Y,tab);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(strcat('Contour Plot: Solution by five-point stencil for ', str));
    saveas(gcf,sprintf('q1_af%d.png',fig));
    fig = fig + 1;
    figure;
    surf(X,Y,error);
    xlabel('x');
    ylabel('y');
    zlabel('Error');
    title(strcat('Absolute Error: five-point stencil for ', str));
    saveas(gcf,sprintf('q1_af%d.png',fig));
    fig = fig + 1;
end

function max_error = five_point_stencil_matrix(h, xd, yd, nx, ny, f, ex, eqx1, eqx2, eqy1, eqy2)
    A = zeros(nx*ny, nx*ny);
    B = zeros(nx*ny, 1);
    for i=1:ny   %at x=0%
        A(i, i)=1;  %u11, u12, u13%
        B(i) = eqy1(yd(i));  
    end

    for i=1:ny  %at x=1%
        A((nx-1)*ny+i,(nx-1)*ny+i) = 1; %u41 u42%
        B((nx-1)*ny+i) = eqy2(yd(i));
    end

    for i=1:nx  %at y=0%
        A((i-1)*ny+1, (i-1)*ny+1) = 1;%u11 u21 u31%
        B((i-1)*ny+1) = eqx1(xd(i));
    end

    for i=1:nx %at y=4%
        A((i)*ny, (i)*ny) = 1;%u14 u24 %
        B((i)*ny) = eqx2(xd(i));
    end

    for i=2:nx-1
        for j=2:ny-1
            A((i-1)*ny+j, (i-1)*ny+j) = -4;
            A((i-1)*ny+j, (i-1)*ny+j-1) = 1;
            A((i-1)*ny+j, (i-1)*ny+j+1) = 1;
            A((i-1)*ny+j, (i)*ny+j) = 1;
            A((i-1)*ny+j, (i-2)*ny+j) = 1;
            B((i-1)*ny+j) = h*h*f(xd(i), yd(j));
        end
    end

    u = A\B;
    tab = zeros(nx, ny);
    tabex = zeros(nx, ny);
    error = zeros(nx, ny);
    max_error=0;
    for i=1:nx
        for j=1:ny
            tab(i, j) = u((i-1)*ny+j);
            tabex(i, j) = ex(xd(i), yd(j));
            error(i, j) = abs(tab(i, j)-tabex(i, j));
            max_error = max([max_error error(i, j)]);
        end
    end
end



