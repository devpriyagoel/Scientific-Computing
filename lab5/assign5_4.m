clear all
clc
clf
close all

obj.t = 0.3:0.02:0.44;
obj.y = sin(obj.t);
obj.n = size(obj.t, 2);
obj.h = obj.t(2:obj.n)-obj.t(1:obj.n-1);
obj.b = (obj.y(2:obj.n)-obj.y(1:obj.n-1))./(obj.h);
obj.v = 2*(obj.h(2:obj.n-1)+obj.h(1:obj.n-2));
obj.u = 6*(obj.b(2:obj.n-1)-obj.b(1:obj.n-2));
obj.z = zeros(obj.n, 1);
obj.u = obj.u';
obj.M = zeros(obj.n-2, obj.n-2);

for i=1:obj.n-2
    obj.M(i, i) = obj.v(i);
end
for i=2:obj.n-2
    obj.M(i-1, i) = obj.h(i);
    obj.M(i, i-1) = obj.h(i);
end
obj.z(2:obj.n-1) = obj.M\obj.u;
obj.S = cell(obj.n-1, 1);
syms x
for i=1:obj.n-1
    obj.S{i} = matlabFunction((obj.z(i+1)*((x-obj.t(i))^3))/(6*obj.h(i))+(obj.z(i)*((obj.t(i+1)-x)^3))/(6*obj.h(i)) +(obj.y(i+1)/obj.h(i)-(obj.z(i+1)*obj.h(i))/6)*(x-obj.t(i)) +(obj.y(i)/obj.h(i)-(obj.z(i)*obj.h(i))/6)*(obj.t(i+1)-x));
end

query.qx = 0.3102;
query.qi = floor((query.qx-obj.t(1))/(obj.t(2)-obj.t(1)) +1);
query.n = size(query.qx, 2);
query.ans = zeros(1, query.n);
query.actual = sin(query.qx);
for i=1:query.n
    query.ans(i) = obj.S{query.qi(i)}(query.qx(i));
end
disp('Value of the interpolant at x = 0.3102 is : ');
disp(query.ans);
disp('Absolute error : ')
disp(abs(query.actual-query.qx));
graph.qx = 0.3:0.001:0.44;
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
plot(graph.qx, sin(graph.qx), 'Linewidth', 3, 'color', 'k');
plot(graph.qx, graph.ans, 'Linewidth', 1.5, 'color', 'y');
title('Natural Spline Interpolant');
legend('Actual function', 'Interpolant');
hold off;
