clear all;
clc;

h=0.05;
k=0.01;
a=0;
b=1;
xd=a:h:b;
td=a:k:b;
nx=(b-a)/h+1; 
nt=(b-a)/k+1;
lambda=2*k/h;
u=zeros(nx,nt);
%icf
icf=@(x) 1 + sin(2*pi*x);
bc_xa=@(x)1;
bc_xb=@(x)1;
u(:,1)=icf(xd);
u(1,:)=bc_xa(td);
u(nx,:)=bc_xb(td);
fig =1;
for j=2:nt
    for i=2:nx-1
            u(i,j)=u(i,j-1) -lambda*(u(i+1,j-1)-u(i-1,j-1))/2 + lambda*lambda*(u(i+1,j-1) - 2*u(i,j-1) + u(i-1,j-1))/2; ;
    end
end
[X,Y]=meshgrid(td,xd);
figure
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title('Approximate Solution by Lax Wendroff (i)');
saveas(gcf,sprintf('q4_1f%d.png',fig));
fig = fig+1;

%b

h=0.05;
k=0.01;
a=0;
b=1;
xd=[a:h:b];
td=[a:k:b];
nx=(b-a)/h+1; 
nt=(b-a)/k+1;
lambda=k/h;
a=2;
u=zeros(nx,nt);
%icf
icf=@(x) 1 + sin(2*pi*x);
bc_xa=@(x)1;
bc_xb=@(x)1;
u(:,1)=icf(xd);
u(1,1)=1;
u(nx,:)=bc_xb(td);

for j=2:nt
    for i=2:nx-1
        u(i,j)=u(i,j-1) - lambda*(u(i+1,j-1)-u(i-1,j-1))/2 + lambda*lambda*(u(i+1,j-1) - 2*u(i,j-1) + u(i-1,j-1))/2; 
    end
    u(1,j)=u(2,j);
end
[X,Y]=meshgrid(td,xd);
figure
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title('Approximate Solution by Lax Wendroff (ii)');
saveas(gcf,sprintf('q4_1f%d.png',fig));
fig = fig+1;
%c

h=0.05;
k=0.01;
a=0;
b=1;
xd=[a:h:b];
td=[a:k:b];
nx=(b-a)/h+1;
nt=(b-a)/k+1;
lambda=k/h;
a=2;
u=zeros(nx,nt);
%icf
icf=@(x) 1 + sin(2*pi*x);
bc_xa=@(x)1;
bc_xb=@(x)1;
u(:,1)=icf(xd);
u(1,1)=1;
u(nx,:)=bc_xb(td);

for j=2:nt
    for i=2:nx-1
        u(i,j)=u(i,j-1) - lambda*(u(i+1,j-1)-u(i-1,j-1))/2 + lambda*lambda*(u(i+1,j-1) - 2*u(i,j-1) + u(i-1,j-1))/2;
    end
    u(1,j)=u(1,j-1) - lambda*(u(2,j-1) - u(1,j-1));
end
[X,Y]=meshgrid(td,xd);
figure
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title('Approximate Solution by Lax Wendroff (iii)');
saveas(gcf,sprintf('q4_1f%d.png',fig));
fig = fig+1;

%d

h=0.05;
k=0.01;
a=0;
b=1;
xd=[a:h:b];
td=[a:k:b];
nx=(b-a)/h+1; 
nt=(b-a)/k+1;
lambda=k/h;
a=2;
u=zeros(nx,nt);
%icf
icf=@(x) 1 + sin(2*pi*x);
bc_xa=@(x)1;
bc_xb=@(x)1;
u(:,1)=icf(xd);
u(1,1)=1;
u(nx,:)=bc_xb(td);

for j=2:nt
    for i=2:nx-1
        u(i,j)=u(i,j-1) - lambda*(u(i+1,j-1)-u(i-1,j-1))/2 + lambda*lambda*(u(i+1,j-1) - 2*u(i,j-1) + u(i-1,j-1))/2; ;
    end
    
    u(1,j)=(u(1,j-1) + u(2,j-1) + lambda*(u(2,j-1) - u(1,j-1)) - u(2,j)*(1+lambda))/(1-lambda);
end
[X,Y]=meshgrid(td,xd);
figure
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title('Approximate Solution by Lax Wendroff (iv)');
saveas(gcf,sprintf('q4_1f%d.png',fig));