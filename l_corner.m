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


%fit data with smoothing splines
% pprho = fit(lambda,lrho,fitType,fitOptions);
if ~isempty(p)
    ppx = csaps(lambda,lrho,p);
    % ppeta = fit(lambda,leta,fitType,fitOptions);
    ppy = csaps(lambda,leta,p);
    %pp-form, first and second derivatives
else
    ppx = csaps(lambda,lrho);
    % ppeta = fit(lambda,leta,fitType,fitOptions);
    ppy = csaps(lambda,leta);
    %pp-form, first and second derivatives
end
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
    figure(fig+1);loglog(rho,eta,rho_c,eta_c,'o');
end

end
