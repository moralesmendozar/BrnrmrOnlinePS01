% for eta
% NB these have not done!
%ODEs
clear; clc;
%odeexamples
options = optimset('Display', 'off');
init = 0.01;
etas = (init:0.01:0.99)'; %this is eta
n = length(etas) ;
y = nan(n);

% define konstants
kappa = 10;
a = .11;
rho = .05;
sigma = .1;
delta = 0.001; % try this, change so get the same rf
rho_bar = 0.02;
q =@(eta) (a + 1/kappa)./(rho.*eta + rho_bar.*(1-eta) + 1/kappa);
iota =@(eta) (q(eta)-1)./kappa;
theta =@(eta) 1 - 1./eta;
%get sigma q2 analytically. Test with numerical too?
sigma_q =@(eta) -(1-eta).*kappa.*(rho - rho_bar)/(1+kappa.*rho).*sigma;
%zeta = (1 - theta).*(sigma + sigma_q);
mu_eta_abs =@(eta) (a - iota(eta))./q(eta) - rho + theta(eta).^2.*(sigma + sigma_q(eta)).^2 .* eta;
sigma_eta_abs =@(eta) -theta(eta).*(sigma + sigma_q(eta)).*eta;
% this next is quite messy, but got by ito on eta so should be fine.
mu_q =@(eta) -(rho - rho_bar).*(a+ 1/kappa)...
    ./(eta.*(rho-rho_bar) + rho_bar + 1/kappa).^2 ...
.*mu_eta_abs(eta)  + kappa^2.*(rho - rho_bar).^2 ...
./(kappa.*(rho.*eta + rho_bar.*(1-eta)) + 1).^2 .* sigma_eta_abs(eta).^2;

% Explicit Euler method to get y
% Works because it satisfies the ODE, and include the boundary constraint,
% so it knows where to start.

% y is D now. D is what we solve for! Unknown now. Will be pdf*sigma^2

xx = etas;
y_explicitEuler = ones(n,1)*init;
y_implicitEuler = ones(n,1)*init;
g=@(eta,y) 2.*mu_eta_abs(eta)./(sigma_eta_abs(eta)).^2.*y;   %this is y'
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
ode = diff(z,t) == 2.*mu_eta_abs(t)./(sigma_eta_abs(t)).^2.*z;
cond = z(init) == init;
ySol(t) = dsolve(ode,cond);

figure(21)
plot(xx,y_explicitEuler,'--r','Linewidth',3)
xlabel('eta')
ylabel('D')
hold on;
plot(xx,y_implicitEuler,'-*b')
plot(xx,ySol(xx),'--g','Linewidth',2)
legend('Explicit Euler','Implicit Euler','Matlab Solution');
hold off;

%%
% BC's
y(1) = init;
g = nan(size(etas));
for i=1:size(etas)-1
    g(i) = 2*mu_eta_abs(i)/sigma_eta_abs(i)^2*y(i); %ODE
    y(i+1) = y(i) + g(i)*(eta(i+1)-eta(i));
end
plot(etas,y)
xlabel('\eta')
ylabel('D')
pdf_notnormalised = y./sigma_eta_abs.^2;
%normalise
pdf = pdf_notnormalised./sum(pdf_notnormalised);
plot(etas,pdf)
xlabel('\eta')
ylabel('pdf')
axis([0 1 0 inf])
saveas(gcf,'p2part4a.pdf')
% save results for y
pdfex = pdf;

%My Explicit Euler method to get pdf(q)
% Works because it satisfies the ODE, and include the boundary constraint,
% so it knows where to start.

% y is D now. D is what we solve for! Unknown now. Will be pdf*sigma^2
y = nan(size(eta));
% BC's
y(1) = 0.01;
for i=1:size(q)-1
    g(i) = 2*mu_q(i)*q(i)/sigma_q(i)^2/q(i)^2*y(i);%ODE
    y(i+1) = y(i) + g(i)*(q(i+1)-q(i));
end



figure(1)
plot(q,y)
xlabel('q')
ylabel('D')
pdf_notnormalised = y./sigma_q.^2./q.^2;
%normalise
pdf = pdf_notnormalised./sum(pdf_notnormalised);

figure(2)
plot(q,pdf)
xlabel('q')
ylabel('pdf')
axis([-inf inf 0 inf])
saveas(gcf,'p2part4a.pdf')
% save results for y
pdfex = pdf;


%My implicit Euler method to get y
% ie need to feed in a proposed solution ie a vector of y's and get a loss
% value out.

%naive solution guess
y = ones(size(eta)).*0.001;

%only optimise on y(2:end)
ysolve = y(2:end);
options = optimoptions('fsolve','Display','iter');
[yimp,fval] = fsolve(@myFunEta,ysolve,options);
yimp = [0; yimp];

figure(3)
plot(eta,yimp)
xlabel('\eta')
ylabel('D')
pdf_notnormalised = yimp./sigma_eta_abs.^2;
%normalise
pdf = pdf_notnormalised./sum(pdf_notnormalised);

figure(4)
plot(eta,pdf)
xlabel('\eta')
ylabel('pdf')
axis([0 1 0 inf])
% save results for y
pdfim = pdf;
%saveas(gcf,'p2part4b.pdf')

%My implicit Euler method to get q
% ie need to feed in a proposed solution ie a vector of y's and get a loss
% value out.
%naive solution guess
tmp = (0.01:0.001:0.99)';
y = ones(size(tmp)).*0.001;
%only optimise on y(2:end)
ysolve = y(2:end);
options = optimoptions('fsolve','Display','iter');
[yimp,fval] = fsolve(@myFunQ,ysolve,options);
yimp = [0; yimp];

figure(5)
plot(q,yimp)
xlabel('\eta')
ylabel('D')
pdf_notnormalised = yimp./sigma_eta_abs.^2;
%normalise
pdf = pdf_notnormalised./sum(pdf_notnormalised);

figure(6)
plot(q,pdf)
xlabel('q')
ylabel('pdf')
axis([0 1 0 inf])
% save results for y
pdfim = pdf;
%saveas(gcf,'p2part4b.pdf')


%use matlab to solve
y0 = 0.01;
[xmatlab,ymatlab] = ode45(@(eta,y)[2.*mu_eta_abs./sigma_eta_abs.^2.*y], eta, y0);
plot(eta,ymatlab)
xlabel('\eta')
ylabel('D')
pdf_notnormalised = ymatlab./sigma_eta_abs.^2;
%normalise
pdf = pdf_notnormalised./sum(pdf_notnormalised);

figure(7)
plot(eta,pdf)
xlabel('\eta')
ylabel('pdf')
axis([0 1 0 inf])
% save results for y
pdfmat = pdf;

%do analytical
yanal = nan(101,1);
yanal(1) = 1;
for i=2:size(x)
    yanal(i) = cos(x(i));
end

figure(8)
plot(x,ymatlab,x,yanal)
legend('matlab','exact')