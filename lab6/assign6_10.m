clc;
clear all;
close all;

length_of_ellipse_arc = @(a,b) (pi.*(2.^(0.5))*((a.*a+b.*b).^(0.5)));
ellipse = @(x,y) ((x.^2)/(9))+(y.^2)/(4);

f = @(x) 4.*((((3.*cos(x)).^2)+((2.*sin(x)).^2)).^(0.5));

len = simpson(0,pi/2,f);

fprintf('Approximate value of the of the ellipse arc’s length evaluated using Simpson rule: %0.4f.\n',len);

function len = simpson(a,b,f)
    n = 1;
    actual_integral = integral(f,a,b);
    len = (b-a)/6*(f(a)+4*f((a+b)/2)+f(b));
    while (abs(real((actual_integral-len)))>=1e-5)
        n = n+1;
        h = (b-a)/n;
        partition = a+h/2:h:b-h/2;
        fx = sum(f(partition-h/2) + 4.*f(partition) + f(partition+h/2));
        len = h/6*fx;
    end
end