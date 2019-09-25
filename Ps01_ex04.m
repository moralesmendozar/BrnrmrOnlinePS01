% UPENN
% Brunnermeier online Course / Princeton
% September 22, 2019.
% Problem set 01, ex 04
%close all;
clear;
clc;
n = 500;
options = optimset('Display', 'off');


%% ex 04_c
lambda = -10;
deltaxs = [0.05 0.1];
%deltaxs = 0.05;
xx = 0:0.05:10;
g=@(x,y) lambda*y;   %this is y'
yreal = exp(lambda.*xx);  %analytical solution
%use matlab to solve the diff equation
syms z(t);
ode = diff(z,t) == lambda*z;
cond = z(0) == 1;
ySol(t) = dsolve(ode,cond);

figure(11)
plot(xx,ySol(xx),'--g','Linewidth',2)
xlabel('x')
ylabel('y')
hold on;
plot(xx,yreal,'-g')

for deltax = deltaxs
    xx = 0:deltax:10;
    n = length(xx);
    y_explicitEuler = ones(n,1);
    y_implicitEuler = ones(n,1);
    for counter = 2:n
        xi = xx(counter-1);
        yi  = y_explicitEuler(counter-1);
        y_explicitEuler(counter) = yi + g(xi,yi)*deltax;
        f=@(y) (y-yi)/(deltax) - g(xx(counter),y);
        [ysolve,fval] = fsolve(@(y)f(y),y_explicitEuler(counter),options);
        y_implicitEuler(counter) = ysolve;
        %y_implicitEuler(counter) = yi/(1-lambda*deltax);
    end
    plot(xx,y_explicitEuler,'--k','Linewidth',1);
    plot(xx,y_implicitEuler,'-b')
end
title('Solution for \Delta x = [0.05 .1]')
legend('Matlab Solution','Analytical solution','Explicit Euler','Implicit Euler');
hold off;


%% ex 04_c
%unstable side
lambda = -10;
deltaxs = [0.19 0.2];
%deltaxs = 0.05;
xx = 0:0.05:10;
g=@(x,y) lambda*y;   %this is y'
yreal = exp(lambda.*xx);  %analytical solution
%use matlab to solve the diff equation
syms z(t);
ode = diff(z,t) == lambda*z;
cond = z(0) == 1;
ySol(t) = dsolve(ode,cond);

figure(12)
plot(xx,ySol(xx),'--g','Linewidth',2)
xlabel('x')
ylabel('y')
hold on;
plot(xx,yreal,'-g')

for deltax = deltaxs
    xx = 0:deltax:10;
    n = length(xx);
    y_explicitEuler = ones(n,1);
    y_implicitEuler = ones(n,1);
    for counter = 2:n
        xi = xx(counter-1);
        yi  = y_explicitEuler(counter-1);
        y_explicitEuler(counter) = yi + g(xi,yi)*deltax;
        f=@(y) (y-yi)/(deltax) - g(xx(counter),y);
        [ysolve,fval] = fsolve(@(y)f(y),y_explicitEuler(counter),options);
        y_implicitEuler(counter) = ysolve;
        %y_implicitEuler(counter) = yi/(1-lambda*deltax);
    end
    plot(xx,y_explicitEuler,'--k','Linewidth',1);
    plot(xx,y_implicitEuler,'-b')
end
title('Solution for \Delta x = [0.19 .2]')
legend('Matlab Solution','Analytical solution','Explicit Euler','Implicit Euler');
hold off;

%% ex 04_c
%unstable side
lambda = -10;
deltaxs = [.21 .25];
%deltaxs = 0.05;
xx = 0:0.05:10;
g=@(x,y) lambda*y;   %this is y'
yreal = exp(lambda.*xx);  %analytical solution
%use matlab to solve the diff equation
syms z(t);
ode = diff(z,t) == lambda*z;
cond = z(0) == 1;
ySol(t) = dsolve(ode,cond);

figure(13)
plot(xx,ySol(xx),'--g','Linewidth',2)
xlabel('x')
ylabel('y')
hold on;
plot(xx,yreal,'-g')

for deltax = deltaxs
    xx = 0:deltax:10;
    n = length(xx);
    y_explicitEuler = ones(n,1);
    y_implicitEuler = ones(n,1);
    for counter = 2:n
        xi = xx(counter-1);
        yi  = y_explicitEuler(counter-1);
        y_explicitEuler(counter) = yi + g(xi,yi)*deltax;
        f=@(y) (y-yi)/(deltax) - g(xx(counter),y);
        [ysolve,fval] = fsolve(@(y)f(y),y_explicitEuler(counter),options);
        y_implicitEuler(counter) = ysolve;
        %y_implicitEuler(counter) = yi/(1-lambda*deltax);
    end
    plot(xx,y_explicitEuler,'--k','Linewidth',1);
    plot(xx,y_implicitEuler,'-b')
end
title('Solution for \Delta x = [0.21 .25]')
legend('Matlab Solution','Analytical solution','Explicit Euler','Implicit Euler');
hold off;
