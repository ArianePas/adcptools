function [intersect] = expon(lambda,rho,p)
    
rho = rho(:);
lambda = lambda(:);

assert(length(rho)==length(lambda),'rho and lambda should be the same length.');



%fit data with smoothing splines
% pprho = fit(lambda,lrho,fitType,fitOptions);
if ~isempty(p)
    ppx = csaps(lambda,rho,p);
    
else
    ppx = csaps(lambda,rho);
    
end
% assume lambda and rho_raw exist



dppx = fnder(ppx);
ddppx = fnder(ppx,2);

line = fnval(ppx,0)+fnval(dppx,0)*lambda;
intersect = -fnval(ppx,0)/fnval(dppx,0);

figure(308)
fnplt(ppx)
hold on
plot(lambda,rho,'.')
plot(lambda(1:5),line(1:5))
plot(intersect,0,'o')
hold off


end