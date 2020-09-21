clc;
clear all;
close all;

syms x;
Q1 = solve_all(((2.*x)./(x.*x-4)), 1, 1.6);
Q2 = solve_all((exp(3.*x).*sin(2.*x)), 0, pi/4 );
Q3 = solve_all(((sin(x)).^2-2.*sin(x).*x+1), 0.75, 1.3);
Q4 = solve_all( 1./(x.*log(x)), exp(1), exp(1)+1);

function Q = solve_all(ff, a, b)
    syms x;
    Q.f =matlabFunction(ff);
    Q.f_dash = matlabFunction(diff(ff));
    Q.f_double_dash = matlabFunction(diff(ff, 2));
    Q.f_triple_dash = matlabFunction(diff(ff, 3));
    Q.f_quad_dash = matlabFunction(diff(ff, 4));
    
    Q.a =a;
    Q.b =b;
    Q.actual = integral(Q.f, Q.a, Q.b);


    Q.rect_left = Q.f(Q.a)*(Q.b-Q.a);
    Q.rect_mid = Q.f((Q.a+Q.b)/2)*(Q.b-Q.a);
    Q.trapezoid = (Q.f(Q.a)+Q.f(Q.b))*(Q.b-Q.a)/2;
    Q.simpsons = ((Q.b-Q.a)/6)*(Q.f(Q.a)+Q.f(Q.b)+4*Q.f(0.5*(Q.a+Q.b)));
    Q.corr_trap = ((Q.b-Q.a)/2)*(Q.f(Q.a)+Q.f(Q.b))-((((Q.b-Q.a).^2)/12).*(Q.f_dash(Q.b)-Q.f_dash(Q.a)));
    Q.simp3by8 = ((Q.b-Q.a)/8)*(Q.f(Q.a)+3*Q.f((Q.b-Q.a)/3)+3*(Q.f(2*(Q.b-Q.a)/3))+Q.f(Q.b));

    Q.rect_actual_error = abs(Q.actual-Q.rect_left);
    Q.rect_mid_actual_error = abs(Q.actual-Q.rect_mid);
    Q.trap_actual_error = abs(Q.actual-Q.trapezoid);
    Q.simp_actual_error = abs(Q.actual-Q.simpsons);
    Q.corr_trap_actual_error = abs(Q.actual-Q.corr_trap);
    Q.simp3by8_actual_error = abs(Q.actual-Q.simp3by8);

    Q.rect_left_error = max(abs(Q.f_dash(Q.a:0.01:Q.b)))*((Q.b-Q.a)^2)/2;
    Q.rect_mid_error = max(abs(Q.f_double_dash(Q.a:0.01:Q.b)))*((Q.b-Q.a)^3)/24;
    Q.trap_error = max(abs(Q.f_double_dash(Q.a:0.01:Q.b)))*((Q.b-Q.a)^3)/12;
    Q.simpsons_error =  max(abs(Q.f_quad_dash(Q.a:0.01:Q.b)))*((Q.b-Q.a)^5)/2880;
    Q.corr_trap_error = max(abs(Q.f_triple_dash(Q.a:0.01:Q.b)))*((((Q.b-Q.a)).^4)/(180));
    Q.simp3by8_error = max(abs(Q.f_quad_dash(Q.a:0.01:Q.b)))*((3.*((((Q.b-Q.a)/3)).^5))/(80));

    for i=1:1
        disp('Numerical Integrals')
        disp('Rectangle Method');
        fprintf('Value: %d, Error Bound: %d, Actual Error: %d\n', Q.rect_left, Q.rect_left_error, Q.rect_actual_error);
        fprintf('Difference in error bound and actual error: %d\n', abs(Q.rect_left_error-Q.rect_actual_error));
        disp('Rectangle Mid-point Method');
        fprintf('Value: %d, Error Bound: %d, Actual Error: %d\n', Q.rect_mid, Q.rect_mid_error, Q.rect_mid_actual_error);
        fprintf('Difference in error bound and actual error: %d\n', abs(Q.rect_mid_error-Q.rect_mid_actual_error));
        disp('Trapezoidal Method');
        fprintf('Value: %d, Error Bound: %d, Actual Error: %d\n', Q.trapezoid, Q.trap_error, Q.trap_actual_error);
        fprintf('Difference in error bound and actual error: %d\n', abs(Q.trap_error-Q.trap_actual_error));
        disp('Corrected Trapezoid Method');
        fprintf('Value: %d, Error Bound: %d, Actual Error: %d\n', Q.corr_trap, Q.corr_trap_error, Q.corr_trap_actual_error);
        fprintf('Difference in error bound and actual error: %d\n', abs(Q.corr_trap_error-Q.corr_trap_actual_error));
        disp('Simpsons Method');
        fprintf('Value: %d, Error Bound: %d, Actual Error: %d\n', Q.simpsons, Q.simpsons_error, Q.simp_actual_error);
        fprintf('Difference in error bound and actual error: %d\n', abs(Q.simpsons_error-Q.simp_actual_error));
        disp('Simpson’s Three-Eighth Method');
        fprintf('Value: %d, Error Bound: %d, Actual Error: %d\n', Q.simp3by8, Q.simp3by8_error, Q.simp3by8_actual_error);
        fprintf('Difference in error bound and actual error: %d\n', abs(Q.simp3by8_error-Q.simp3by8_actual_error));
    end
end
                                                                           