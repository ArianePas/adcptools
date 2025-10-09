function animate_solution(F, X, name, sav)
fi = figure;
filename = name;
ncolor = 100;
amax = max(abs([min(F, [], "all"), max(F, [], "all")]));
levels = linspace(-amax, amax, ncolor);
[~,ha]=contourf(squeeze(X.Y(1,:,:))', squeeze(X.Z(1,:,:))',...
    squeeze(F(1,:,:))' , levels, "LineColor",'none');
colorbar;
title(name)
cm=colormap(gca, helpers.cmaps("velmap"));
clim([-amax, amax])
ylim([min(X.Z, [], 'all'), max(X.Z, [], 'all')])
set(gca, 'XDir','reverse') % Very important
for tim = 1:1:size(F,1)
    frame = getframe(fi);
    im = frame2im(frame);
    ha.YData = squeeze(X.Z(tim,:,:))';
    ha.ZData = squeeze(F(tim,:,:))';
    if tim == 1
        [imind,cm] = rgb2ind(im,256);
        if sav
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
        end
    else
        imind = rgb2ind(im, cm);
        if sav
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end
    end
    pause(.2)
end
end