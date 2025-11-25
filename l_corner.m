function [rho_c,eta_c,lambda_c] = l_corner(rho,eta,lambda,p,fig)
% Finds the corner of the rho and eta L-curve: trade-off between the size 
% of the regularized solutions and their fit to the data, as the regularization 
% parameter varies, displyed on a log-log plot. 
% 
% Parameters
% ==========
%    rho: Lx1 vector, list of ridge residual norm (i.e. rmse) computed for  a set of lambda values 
%    eta: Lx1 vector, list of norm of regularized ridge solutions (i.e. beta-beta_prior distance) 
%         computed for a set of lambda values 
% lambda: Lx1 vector of lambda values
%      p: number, cubic spline smoothness for interpolating the l-curve
%    fig: int or [], figure number for curve plotting. If empty, don't plot the results
% 
% Returns
% =======
%    rho_c: float, rho value with maximum curvature (i.e. for the best lambda)
%    eta_c: float, eta value with maximum curvature (i.e. for the best lambda)
% lambda_c: float, best lambda (i.e. value with maximum curvature)
% 
%     
    
rho = rho(:);
eta = eta(:);
lambda = lambda(:);

assert(length(rho)==length(lambda),'rho and lambda should be the same length.');
assert(length(eta)==length(lambda),'eta and lambda should be the same length.');

lrho = log(rho);
leta = log(eta);

p = 0.9


%fit data with smoothing splines
% pprho = fit(lambda,lrho,fitType,fitOptions);
if ~isempty(p)
    ppx = csaps(lambda,lrho,p); %fit smoothing spline
    
    % ppeta = fit(lambda,leta,fitType,fitOptions);
    ppy = csaps(lambda,leta,p);
    %pp-form, first and second derivatives
else
    ppx = csaps(lambda,lrho);
    % ppeta = fit(lambda,leta,fitType,fitOptions);
    ppy = csaps(lambda,leta);
    %pp-form, first and second derivatives
end
% assume lambda and rho_raw exist
% 1) clean data and compute lrho (natural log)
valid = isfinite(lambda) & isfinite(rho_raw) & (rho_raw > 0);
lambda = lambda(valid);
rho_raw = rho_raw(valid);
lrho = log(rho_raw);    % natural log

% 2) unique/sort lambda
[lambda, idx] = unique(lambda, 'stable');   % or 'sorted'
lrho = lrho(idx);

% 3) fit csaps in log-domain
p = 0.9;
pp = csaps(lambda, lrho, p);

% 4) evaluate on fine grid and back-transform
lam_grid = linspace(min(lambda), max(lambda), 1000);
lrho_fit = fnval(pp, lam_grid);
rho_fit = exp(lrho_fit);    % back-transform

% 5) optional bias correction using residual variance
lrho_at_data = fnval(pp, lambda);
res = lrho - lrho_at_data;
sigma2 = mean(res.^2);      % unbiased: use var(res,1) or var(res,0) as needed
rho_fit_biascorr = exp(lrho_fit + 0.5*sigma2);

% 6) plots
figure;
subplot(2,1,1);
plot(lambda, lrho, 'ko'); hold on;
plot(lam_grid, lrho_fit, 'b-','LineWidth',1.2);
xlabel('\lambda'); ylabel('ln(\rho)'); legend('data','csaps fit');

subplot(2,1,2);
plot(lambda, rho_raw, 'ko'); hold on;
plot(lam_grid, rho_fit, 'b-','LineWidth',1.2);
plot(lam_grid, rho_fit_biascorr, 'r--','LineWidth',1);
set(gca,'YScale','log'); % optional: log-y for clarity
xlabel('\lambda'); ylabel('\rho'); legend('data','exp(fit)','bias-corr');
dppx = fnder(ppx);
dppy = fnder(ppy);
ddppx = fnder(ppx,2);
ddppy = fnder(ppy,2);

%evaluate functions between lambda min and max
dfrho = fnval(dppx,lambda);
dfeta = fnval(dppy,lambda);
ddfrho = fnval(ddppx,lambda);
ddfeta = fnval(ddppy,lambda);

%cuvature
kappa = (dfrho.*ddfeta - ddfrho.*dfeta)./(dfrho.^2 + dfeta.^2).^1.5;

%maximum curvature
% [~,ikappamax] = max(abs(kappa));
[~,ikappamax] = max(kappa);

%---TEST---%
% zrho = fnzeros(dppx);
% if ~isempty(zrho)
%     lambda_c = zrho(1,end);
% else
%     lambda_c = lambda(1);
% end
% rho_c = exp(fnval(ppx,lambda_c));
% eta_c = exp(fnval(ppy,lambda_c));
% kappa_c = interp1(lambda,kappa,lambda_c,'spline');
% if nargin>4
%     figure(fig);semilogx(lambda,kappa,lambda_c,kappa_c,'o');
%     figure(fig+1);loglog(rho,eta,rho_c,eta_c,'o');
% end
% return;
%----------%

rho_c = rho(ikappamax);
eta_c = eta(ikappamax);
lambda_c = lambda(ikappamax);

% ipt = findchangepts(leta,'Statistic','linear','MaxNumChanges',2);
% iptm = round(mean(ipt));
% rho_c = rho(iptm);
% eta_c = eta(iptm);
% lambda_c = lambda(iptm);

if nargin>4
    figure(fig);semilogx(lambda,kappa,lambda_c,kappa(ikappamax),'o');
    figure(fig+1);loglog(rho,eta, '.' ,rho_c,eta_c,'o');
end

end
