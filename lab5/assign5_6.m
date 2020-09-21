clear all
clc
clf
close all

obj.t = [0 3 5 8 13];
obj.y = [0 225 383 623 993];
obj.n = size(obj.t, 2);
obj.h = obj.t(2:obj.n)-obj.t(1:obj.n-1);
nat.b = (obj.y(2:obj.n)-obj.y(1:obj.n-1))./(obj.h);
nat.v = 2*(obj.h(2:obj.n-1)+obj.h(1:obj.n-2));
nat.u = 6*(nat.b(2:obj.n-1)-nat.b(1:obj.n-2));
nat.z = zeros(obj.n, 1);
nat.u = nat.u';
nat.M = zeros(obj.n-2, obj.n-2);

for i=1:obj.n-2
    nat.M(i, i) = nat.v(i);
end
for i=2:obj.n-2
    nat.M(i-1, i) = obj.h(i);
    nat.M(i, i-1) = obj.h(i);
end
nat.z(2:obj.n-1) = nat.M\nat.u;
nat.S = cell(obj.n-1, 1);
nat.S_der = cell(obj.n-1, 1);
syms x
for i=1:obj.n-1
    nat.S{i} = matlabFunction((nat.z(i+1)*((x-obj.t(i))^3))/(6*obj.h(i))+(nat.z(i)*((obj.t(i+1)-x)^3))/(6*obj.h(i)) +(obj.y(i+1)/obj.h(i)-(nat.z(i+1)*obj.h(i))/6)*(x-obj.t(i)) +(obj.y(i)/obj.h(i)-(nat.z(i)*obj.h(i))/6)*(obj.t(i+1)-x));
    nat.S_der{i} = matlabFunction((nat.z(i+1)*((x-obj.t(i))^2))/(2*obj.h(i))-(nat.z(i)*((obj.t(i+1)-x)^2))/(2*obj.h(i)) +(obj.y(i+1)/obj.h(i)-(nat.z(i+1)*obj.h(i))/6) -(obj.y(i)/obj.h(i)-(nat.z(i)*obj.h(i))/6));
end

clm.b = [75 (obj.y(2:obj.n)-obj.y(1:obj.n-1))./(obj.h) 72];
clm.v = [obj.h(1)/3 (obj.h(2:obj.n-1)+obj.h(1:obj.n-2))/3 obj.h(obj.n-1)/3];
clm.u = clm.b(2:obj.n+1)-clm.b(1:obj.n);
clm.u = clm.u';
clm.M = zeros(obj.n, obj.n);
for i=1:obj.n
    clm.M(i, i) = clm.v(i);
end
for i=2:obj.n
    clm.M(i-1, i) = obj.h(i-1)/6;
    clm.M(i, i-1) = obj.h(i-1)/6;
end

clm.z = clm.M\clm.u;
clm.S = cell(obj.n-1, 1);
clm.S_der = cell(obj.n-1, 1);
syms x
for i=1:obj.n-1
    clm.S{i} = matlabFunction((clm.z(i+1)*((x-obj.t(i))^3))/(6*obj.h(i))+(clm.z(i)*((obj.t(i+1)-x)^3))/(6*obj.h(i)) +(obj.y(i+1)/obj.h(i)-(clm.z(i+1)*obj.h(i))/6)*(x-obj.t(i)) +(obj.y(i)/obj.h(i)-(clm.z(i)*obj.h(i))/6)*(obj.t(i+1)-x));
    clm.S_der{i} = matlabFunction((clm.z(i+1)*((x-obj.t(i))^2))/(2*obj.h(i))-(clm.z(i)*((obj.t(i+1)-x)^2))/(2*obj.h(i)) +(obj.y(i+1)/obj.h(i)-(clm.z(i+1)*obj.h(i))/6) -(obj.y(i)/obj.h(i)-(clm.z(i)*obj.h(i))/6));

end

query.qx = 10;
query.n = size(query.qx, 2);
for i=1:query.n
    query.qi(i) = sum(obj.t<=query.qx(i));
end
query.nat_position = zeros(1, query.n);
query.nat_speed = zeros(1, query.n);
query.clm_position = zeros(1, query.n);
query.clm.speed = zeros(1, query.n);
for i=1:query.n
    if(query.qi(i)==obj.n)
        query.nat_position(i) = 13;
        query.nat_speed(i) = 72;
        query.clm_position(i) = 13;
        query.clm_speed(i) = 72;
        continue;
    end
    query.nat_position(i) = nat.S{query.qi(i)}(query.qx(i));
    query.nat_speed(i) = nat.S_der{query.qi(i)}(query.qx(i));
    query.clm_position(i) = clm.S{query.qi(i)}(query.qx(i));
    query.clm_speed(i) = clm.S_der{query.qi(i)}(query.qx(i));
end

disp('Using Natural Spline Interpolation : ')
fprintf('Position of the car when t = 10 seconds :');
fprintf('%f feets \n', query.nat_position);
fprintf('Speed of the car when t = 10 seconds :');
fprintf('%f feets/sec \n\n', query.nat_speed);
disp('Using Clamped Spline Interpolation : ')
fprintf('Position of the car when t = 10 seconds :');
fprintf('%f feets \n', query.clm_position);
fprintf('Speed of the car when t = 10 seconds :');
fprintf('%f feets/sec \n\n', query.clm_speed);
graph.qx = 0:0.01:13;
graph.n = size(graph.qx, 2);
for i=1:graph.n
    graph.qi(i) = sum(obj.t<=graph.qx(i));
end
graph.nat_position = zeros(1, graph.n);
graph.nat_speed = zeros(1, graph.n);
graph.clm_position = zeros(1, graph.n);
graph.clm.speed = zeros(1, graph.n);
for i=1:graph.n
    if(graph.qi(i)==obj.n)
        graph.nat_position(i) = 993;
        graph.nat_speed(i) = 72;
        graph.clm_position(i) = 993;
        graph.clm_speed(i) = 72;
        continue;
    end
    graph.nat_position(i) = nat.S{graph.qi(i)}(graph.qx(i));
    graph.nat_speed(i) = nat.S_der{graph.qi(i)}(graph.qx(i));
    graph.clm_position(i) = clm.S{graph.qi(i)}(graph.qx(i));
    graph.clm_speed(i) = clm.S_der{graph.qi(i)}(graph.qx(i));
end

figure(1)
hold on;
plot(graph.qx, graph.nat_position);
plot(graph.qx, graph.clm_position);
title('Spline Interpolants for Position of the car');
legend('Natural Spline', 'Clamped Spline');
hold off;
figure(2)
hold on;
plot(graph.qx, graph.nat_speed);
plot(graph.qx, graph.clm_speed);
title('Spline Interpolants for Speed of the car');
legend('Natural Spline', 'Clamped Spline');
hold off;