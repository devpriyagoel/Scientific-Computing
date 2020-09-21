clear all
clc
clf
close all

a = 0;
b = 13.25;

mi_per_hr_to_ft_per_sec =  @(x) x/0.00018939/3600;

time = [0,3,5,8,13]; 
car.pos = [0,225,383,623,993];
car.speed = [75,77,80,74,72];
car.dd_table = divided_difference(time,car.pos,car.speed);
hermite = Hermite_Interpolation(car.dd_table);


% Part (a)

t = 10;
position = hermite.distance(t);
speed = hermite.speed(t);
disp('Hermite Interpolation :');
fprintf('t = %d seconds \nPosition of the car = %f feet\nSpeed of the car = %f feet/second\n\n',t,position,speed);

% Part (b)
limiting_val = 55;
limit = mi_per_hr_to_ft_per_sec(limiting_val);
limit_exceeded = @(x)hermite.speed(x) - limit;
x_not = 5;
limit_exceeded_at = Newton(limit_exceeded,hermite.acceleration, x_not);
fprintf('Speed limit of %d mi/hr is exceeded at %0.14f seconds.\n\n',limiting_val,(limit_exceeded_at));

% Part (c)


x_not = 14;
max_speed_at = Newton(hermite.acceleration,hermite.jerk, x_not);
max_speed = hermite.speed(max_speed_at);
fprintf('Predicted maximum speed of the car is %0.10f ft/s.\n',max_speed);

figure(1);
hold on;
fplot(@(t)hermite.speed(t),[a,b],'color','magenta');
line([a,b],[limit,limit],'color','blue');
line([a,b],[max_speed,max_speed],'color','black');
legend('speed (in ft/s)','speed = 55 mi/hr','maximum speed')
hold off;


function X = divided_difference(x,y,z)
    n = 2*length(x);
    X = zeros(n, n);
    for i = 1:n
        if(mod(i,2)==1)
            X(i,1)=x(ceil(i/2));
            X(i,2)=y(ceil(i/2));
            
        else
            if (i~=2)
                X(i-1,3) = (X(i-1,2)-X(i-2,2))/ (X(i-1,1)-X(i-2,1));
            end
            X(i,1)=x(i/2);
            X(i,2)=y(i/2);
            X(i,3)=z(i/2);
        end
    end
    for j = 4:n+1
        for i = j-1:n
            X(i,j) = (X(i,j-1) - X(i-1,j-1))/(X(i,1) - X(i-j+2,1));
        end
    end
end

function hermite = Hermite_Interpolation(X)
    n = size(X, 1);
    syms x
    for i = 2:n+1
        Y(i-1) = X(1,i);
    end
    hermite.distance = Y(1);
    for i = 1:n-1
        temp = 1;
        for j = 1:i
            temp = temp*(x-X(j,1));
        end
        hermite.distance = hermite.distance+temp*X(i+1,i+2);
    end
    hermite.speed = diff(hermite.distance);
    hermite.acceleration = diff(hermite.speed);
    hermite.jerk = diff(hermite.acceleration);
    hermite.distance = matlabFunction(hermite.distance);
    hermite.speed = matlabFunction(hermite.speed);
    hermite.acceleration = matlabFunction(hermite.acceleration);
    hermite.jerk = matlabFunction(hermite.jerk);
end

function t = Newton(a,a_dash,t)
    t_new = t - a(t)/a_dash(t);
    while(abs((t-t_new))>1e-15)
        t = t_new;
        t_new = t_new - a(t_new)/a_dash(t_new);
    end
end