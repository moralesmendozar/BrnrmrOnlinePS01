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
deltaxs = [0.05 0.1 0.19 0.2 .21 .25];
deltaxs = 0.05;
xx = 0:0.05:10;
yreal = exp(lambda.*xx);  %analytical solution
%use matlab to solve the diff equation
syms z(t);
ode = diff(z,t) == t*cos(t^2)*z^2;
cond = z(0) == 1;
ySol(t) = dsolve(ode,cond);

figure(11)
plot(xx,ySol(xx),'--g','Linewidth',2)
plot(xx,yreal,'-g')
hold on;

for deltax = deltaxs
    xx = 0:deltax:10;
    y_explicitEuler = ones(n,1);
    y_implicitEuler = ones(n,1);
    
    g=@(x,y) lambda*y;   %this is y'
    for counter = 2:n
        xi = xx(counter-1);
        yi  = y_explicitEuler(counter-1);
        y_explicitEuler(counter) = yi + g(xi,yi)*deltax;
        f=@(y) (y-yi)/(deltax) - g(xx(counter),y);
        [ysolve,fval] = fsolve(@(y)f(y),yi,options);
        y_implicitEuler(counter) = ysolve;
    end

    plot(xx,y_explicitEuler,'--k','Linewidth',3);
    plot(xx,y_implicitEuler,'-b')
end

legend('Matlab Solution','Analytical solution','Explicit Euler','Implicit Euler',);
hold off;


%% just sinking in a safe copy hehe

for deltax = deltaxs
    deltax = 0.05
    xx = linspace(0,10,n);
    y_explicitEuler = ones(n,1);
    y_implicitEuler = ones(n,1);
    yreal = 1./(1-sin(xx.^2)./2);
    g=@(x,y) x*cos(x^2)*y^2;   %this is y'
    for counter = 2:n
        xi = xx(counter-1);
        yi  = y_explicitEuler(counter-1);
        y_explicitEuler(counter) = yi + g(xi,yi)*(xx(counter)-xi);
        f=@(y) (y-yi)/(xx(counter)-xi) - g(xx(counter),y);
        [ysolve,fval] = fsolve(@(y)f(y),yi,options);
        y_implicitEuler(counter) = ysolve;
    end
    %use matlab to solve
    syms z(t);
    ode = diff(z,t) == t*cos(t^2)*z^2;
    cond = z(0) == 1;
    ySol(t) = dsolve(ode,cond);

    figure(2)
    plot(xx,y_explicitEuler,'--r','Linewidth',3)
    hold on;
    plot(xx,y_implicitEuler,'-*b')
    plot(xx,ySol(xx),'--g','Linewidth',2)
    plot(xx,yreal,'-k')
    legend('Explicit Euler','Implicit Euler','Matlab Solution','Analytical solution');
    hold off;
end
%% ex 03_c
xx = linspace(0,10,n);
y_explicitEuler = ones(n,1);
y_implicitEuler = ones(n,1);
v = zeros(n,1);  % dy/dx = v,   dv/dt = d^2y/dx^2 = -y
yreal = cos(xx);
g=@(x,y) -y;   %this is y'' = -y

for counter = 2:n
    xi = xx(counter-1);
    yi  = y_explicitEuler(counter-1);
    vi  = v(counter-1);
    v(counter) = vi - yi*(xx(counter)-xi);
    y_explicitEuler(counter) = yi + v(counter)*(xx(counter)-xi);
    %f=@(y) (y-yi)/(xx(counter)-xi) - g(xx(counter),y);
    %[ysolve,fval] = fsolve(@(y)f(y),yi,options);
    %y_implicitEuler(counter) = ysolve;
    y_implicitEuler(counter) = y_explicitEuler(counter);
end
%use matlab to solve
syms z(t);
Dz = diff(z);
ode = diff(z,t,2) == -z;
cond1 = z(0) == 1;
cond2 = Dz(0) == 0;
conds = [cond1 cond2];
ySol(t) = dsolve(ode,conds);

figure(3)
plot(xx,y_explicitEuler,'--r','Linewidth',3)
hold on;
%plot(xx,y_implicitEuler,'-*b')
plot(xx,ySol(xx),'--g','Linewidth',2)
plot(xx,yreal,'-k')
legend('Explicit Euler','Matlab Solution','Analytical solution');
hold off;