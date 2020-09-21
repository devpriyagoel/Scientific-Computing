clc;
clear all;
close all;

part = 0:6:84;
f = [124 134 148 156 147 133 121 109 99 85 78 89 104 116 123];

Simpson = simpson(f,6);
fprintf('Approximate  value of track length evaluated using Simpson rule: %d feet\n',Simpson);

function Simpson = simpson(f,h)
    n = length(f);
    areas = h/3*(f(1:2:n-2)+4*f(2:2:n-1)+f(3:2:n));
    Simpson = sum(areas);
end
