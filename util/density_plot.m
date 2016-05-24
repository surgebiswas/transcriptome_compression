function hc = density_plot(x,y)

    [~,density,X,Y]=kde2d([x, y]);
    rdensity = (density - min(min(density)))/max(max(density));
    [hc, hc] = contourf(X,Y,rdensity,100);
    set(hc, 'LineStyle', 'none')
    colorbar
    caxis([0 1]);



end

