% for eta
% NB these have not done!
%ODEs
clear; clc;
%odeexamples
eta = (0.01:0.01:0.99)'; %this is eta
y = nan(size(eta));
g = nan(size(eta));

% BC's
y(1) = 0.01;

% define konstants
kappa = 10;
a = .11;
rho = .05;
sigma = .1;
delta = 0.0; % try this, change so get the same rf
rho_bar = 0.02;
q = (a + 1/kappa)./(rho.*eta + rho_bar.*(1-eta) + 1/kappa);
iota = (q-1)./kappa;
theta = 1 - 1./eta;
%get sigma q2 analytically. Test with numerical too?
sigma_q = -(1-eta).*kappa.*(rho - rho_bar)/(1+kappa.*rho).*sigma;
%zeta = (1 - theta).*(sigma + sigma_q);
mu_eta_abs = (a - iota)./q - rho + theta.^2.*(sigma + sigma_q).^2 .* eta;
sigma_eta_abs = -theta.*(sigma + sigma_q).*eta;
% this next is quite messy, but got by ito on eta so should be fine.
mu_q = -(rho - rho_bar).*(a+ 1/kappa)...
    ./(eta.*(rho-rho_bar) + rho_bar + 1/kappa).^2 ...
.*mu_eta_abs.*eta.*kappa  + kappa^2.*(rho - rho_bar).^2 ...
./(kappa.*(rho.*eta + rho_bar.*(1-eta)) + 1).^2 .* sigma_eta_abs.^2;

%My Explicit Euler method to get y
% Works because it satisfies the ODE, and include the boundary constraint,
% so it knows where to start.

% y is D now. D is what we solve for! Unknown now. Will be pdf*sigma^2

for i=1:size(eta)-1
    g(i) = 2*mu_eta_abs(i)/sigma_eta_abs(i)^2*y(i); %ODE
    y(i+1) = y(i) + g(i)*(eta(i+1)-eta(i));
end
plot(eta,y)
xlabel('\eta')
ylabel('D')
pdf_notnormalised = y./sigma_eta_abs.^2;
%normalise
pdf = pdf_notnormalised./sum(pdf_notnormalised);
plot(eta,pdf)
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
plot(q,y)
xlabel('q')
ylabel('D')
pdf_notnormalised = y./sigma_q.^2./q.^2;
%normalise
pdf = pdf_notnormalised./sum(pdf_notnormalised);
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
plot(eta,yimp)
xlabel('\eta')
ylabel('D')
pdf_notnormalised = yimp./sigma_eta_abs.^2;
%normalise
pdf = pdf_notnormalised./sum(pdf_notnormalised);
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
plot(q,yimp)
xlabel('\eta')
ylabel('D')
pdf_notnormalised = yimp./sigma_eta_abs.^2;
%normalise
pdf = pdf_notnormalised./sum(pdf_notnormalised);
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


plot(x,ymatlab,x,yanal)
legend('matlab','exact')