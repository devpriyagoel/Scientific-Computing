clear all;
clc;

h=0.005;
k=0.001;
a=0;
b=1;
xd=a:h:b;
td=a:k:b;
nx=(b-a)/h+1;
nt=(b-a)/k+1;
lambda=-2*k/h;
tab=zeros(nx,nt);

icf=@(x) 1 + sin(2*pi*x);
bc_xb=@(x)1;
tab(:,1)=icf(xd);
tab(nx,:)=bc_xb(td);

for i=nx-1:-1:1
    for j=2:nt
    tab(i,j)=(1+lambda)*tab(i,j-1) - lambda*tab(i+1,j-1);
    end
end
fig=1;
[X,Y]=meshgrid(td,xd);
surf(X,Y,tab);
xlabel('x');
ylabel('y');
zlabel('u');
title('Approximate Solution using FTFS');
saveas(gcf,sprintf('q3_1f%d.png',fig));
fig = fig+1;
tab=zeros(nx,nt);
%icf
icf=@(x) 1 + sin(2*pi*x);
bc_xb=@(x)1;
tab(:,1)=icf(xd);
tab(nx,:)=bc_xb(td);
temp=nx-1;
while(temp>0)
    for j=2:nt
    tab(temp,j)=(tab(temp,j-1) - lambda*tab(temp+1,j))*(1/(1-lambda));
    end
    temp=temp-1;
end
[X,Y]=meshgrid(td,xd);
figure
surf(X,Y,tab);
xlabel('x');
ylabel('y');
zlabel('u');
title('Approximate Solution using BTFS');
saveas(gcf,sprintf('q3_1f%d.png',fig));

