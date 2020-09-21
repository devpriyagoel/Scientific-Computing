clear all
clc
clf
close all

obj.t = 0.1:0.1:0.4;
obj.y = [-0.62049958 -0.28398668 0.00660095 0.24842440];
obj.n = size(obj.t, 2);
obj.h = obj.t(2:obj.n)-obj.t(1:obj.n-1);
obj.f = @(x)x.*cos(x)-2*(x.^2)+3*x-1;
obj.b = [3.58502082 (obj.y(2:obj.n)-obj.y(1:obj.n-1))./(obj.h) 2.16529366];
obj.v = [obj.h(1)/3 (obj.h(2:obj.n-1)+obj.h(1:obj.n-2))/3 obj.h(obj.n-1)/3];
obj.u = obj.b(2:obj.n+1)-obj.b(1:obj.n);
obj.z = zeros(obj.n, 1);
obj.u = obj.u';
obj.M = zeros(obj.n, obj.n);
for i=1:obj.n
    obj.M(i, i) = obj.v(i);
end
for i=2:obj.n
    obj.M(i-1, i) = obj.h(i-1)/6;
    obj.M(i, i-1) = obj.h(i-1)/6;
end

obj.z = obj.M\obj.u;
obj.S = cell(obj.n-1, 1);
syms x
for i=1:obj.n-1
    obj.S{i} = matlabFunction((obj.z(i+1)*((x-obj.t(i))^3))/(6*obj.h(i))+(obj.z(i)*((obj.t(i+1)-x)^3))/(6*obj.h(i)) +(obj.y(i+1)/obj.h(i)-(obj.z(i+1)*obj.h(i))/6)*(x-obj.t(i)) +(obj.y(i)/obj.h(i)-(obj.z(i)*obj.h(i))/6)*(obj.t(i+1)-x));
end
query.qx = 0.2013;
query.qi = floor((query.qx-obj.t(1))/(obj.t(2)-obj.t(1)) +1);
query.n = size(query.qx, 2);
query.ans = zeros(1, query.n);


for i=1:query.n
    query.ans(i) = obj.S{query.qi(i)}(query.qx(i));
end
disp('Value of the interpolant at x = 0.2013 is : ');
disp(query.ans);
disp('Error of the interpolant at x = 0.2013 is :');
disp(abs(query.ans-obj.f(query.qx)));

graph.qx = 0.1:0.001:0.4;
graph.n = size(graph.qx, 2);
graph.ans = zeros(1, graph.n);
for i=1:graph.n
    graph.qi(i) = sum(obj.t<=graph.qx(i));
end
for i=1:graph.n
    if(graph.qi(i)==obj.n)
        graph.ans(i)=obj.y(obj.n);
        continue;
    end
    graph.ans(i) = obj.S{graph.qi(i)}(graph.qx(i));
end
figure(1);
hold on;
plot(graph.qx, graph.ans, 'Linewidth', 3, 'color', 'k');
plot(graph.qx, obj.f(graph.qx), 'Linewidth', 1.5, 'color', 'y');
title('Clamped Spline Interpolant');
legend('Actual function', 'Interpolant');
hold off;
