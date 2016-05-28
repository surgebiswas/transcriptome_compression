function [sx, sy, in] = jbm_plot_venn( px, py, cm, plotname )

    if length(px) > 1000 % working with genes. do FDR correction
        st = 0.01;
    else
        st = 0.05;
    end

    figure;
    sx = sum(px < st);
    sy = sum(py < st);
    in = sum(px < st &  py < st);
    h = venn([sx, sy], in);
    set(h(1), 'FaceColor', cm(10,:))
    set(h(2), 'FaceColor', cm(end-10,:))
    axis off
    plotSave(plotname);
    close

end

