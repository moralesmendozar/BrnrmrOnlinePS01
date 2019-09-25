%problem 2.2 brunnermeier class
%cd('C:\Users\anbr.fi\Dropbox (CBS Finance)\Year3\ECO529 Macro, Money and International finance (online course)')
clearvars
% define konstants
kappa = 10;
a = .11;
rho = .05;
sigma = .1;
delta = 0.001; % try this, change so get the same rf
% this only needed for second model
rho_bar = 0.02;

%make eta grid
eta = 0:0.01:1;

%first replicate simple model %nb dynamic vars are t subscribted
q = (1 + kappa*a)/(1 + kappa*rho); %ie constant
r = rho + 1/kappa*log(q) - delta - sigma^2./eta;
%absolute drift eta
mu_eta_abs = (1-eta).^2./eta.^2*sigma^2.*eta;
%absolute vol eta
sigma_eta_abs = (1-eta)./eta*sigma.*eta;

figure(1)
%plot
subplot(2,2,1)
plot(eta,q.*ones(101,1))
axis([0 1 1 2])
xlabel('\eta')
ylabel('q')

subplot(2,2,2)
plot(eta,r)
axis([0 1 -0.4 0.1])
xlabel('\eta')
ylabel('r')

subplot(2,2,3)
plot(eta,sigma_eta_abs)
xlabel('\eta')
ylabel('\eta \sigma^{\eta}')

subplot(2,2,4)
plot(eta,mu_eta_abs)
xlabel('\eta')
ylabel('\eta \mu^{\eta}')

%second model
%clearvars -except kappa a rho sigma delta rho_bar eta
%make new model variables
q2 = (a + 1/kappa)./(rho.*eta + rho_bar.*(1-eta) + 1/kappa);
iota2 = (q2-1)./kappa;
theta2 = 1 - 1./eta;
%get sigma q2 analytically. Test with numerical too?
sigma_q2 = -(1-eta).*kappa.*(rho - rho_bar)/(1+kappa.*rho).*sigma;
zeta2 = (1 - theta2).*(sigma + sigma_q2);
%mu_eta_abs2 = ((a - iota2)./q2 - rho + ((1-eta)./eta).^2.*(sigma + sigma_q2).^2).*eta;
%theta2.^2.*(sigma + sigma_q2).^2 .* eta;
mu_eta_abs2 = ((a - iota2)./q2 - rho + theta2.^2.*(sigma + sigma_q2).^2 ).*eta;
sigma_eta_abs2 = -theta2.*(sigma + sigma_q2).*eta;


% this next is quite messy, but got by ito on eta so should be fine.
%mu_q2 = -(rho - rho_bar).*(a+ 1/kappa)...
%    ./(eta.*(rho-rho_bar) + rho_bar + 1/kappa).^2 ...
%.*mu_eta_abs2.*eta.*kappa  + kappa^2.*(rho - rho_bar).^2 ...
%./(kappa.*(rho.*eta + rho_bar.*(1-eta)) + 1).^2 .* sigma_eta_abs2.^2;
mu_q2 = -(rho - rho_bar).*(a+ 1/kappa)...
    ./(eta.*(rho-rho_bar) + rho_bar + 1/kappa).^2 ...
.*mu_eta_abs2  + kappa^2.*(rho - rho_bar).^2 ...
./(kappa.*(rho.*eta + rho_bar.*(1-eta)) + 1).^2 .* sigma_eta_abs2.^2;
r2 = (a-iota2)./q2 + 1/kappa*log(1+kappa.*iota2) - delta + mu_q2 + ...
    sigma*sigma_q2 - zeta2.*(sigma + sigma_q2);

figure(2)
%plot
subplot(2,2,1)
plot(eta,q.*ones(101,1),eta,q2)
axis([0 1 1 2])
xlabel('\eta')
ylabel('q')
legend('\rho^h = \rho^e','\rho^h < \rho^e')

subplot(2,2,2)
plot(eta,r,eta,r2)
axis([0 1 -0.4 0.1])
xlabel('\eta')
ylabel('r')

subplot(2,2,3)
plot(eta,sigma_eta_abs,eta,sigma_eta_abs2)
xlabel('\eta')
ylabel('\eta \sigma^{\eta}')

subplot(2,2,4)
plot(eta,mu_eta_abs,eta,mu_eta_abs2)
axis([0 1 -0.015 0.1])
xlabel('\eta')
ylabel('\eta \mu^{\eta}')
saveas(gcf,'p2part2.pdf')