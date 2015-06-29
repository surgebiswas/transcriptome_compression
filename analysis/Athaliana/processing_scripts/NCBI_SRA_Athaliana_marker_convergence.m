function NCBI_SRA_Athaliana_marker_convergence( lY, tids, qt, results )
% results struct output from NCBI_SRA_Athaliana_power_analysis.m

NMARKERS = 10;

S = [];
for i = 1 : length(results.somps)
    S = [S; results.somps{i}.S];
end
S = S';

fS = S(:,end);
corrwf = corr(lY, lY(:, fS(1:NMARKERS)))';

corrwfs = corrwf; % Standardize in future?
% for i = 1 : size(corrwfs,1)
%     corrwfs(i,:) = (corrwfs(i,:) - mean(corrwfs(i,:)))./std(corrwfs(i,:));
% end


R = zeros(NMARKERS, size(S,2));
for i = 1 : NMARKERS
    for j = 1 : size(S,2)
        R(i,j) = corrwfs(i, S(i,j));
    end
end

cm = cbrewer('seq', 'YlOrRd', 100, 'cubic');


R(:,end) = []; % remove the last col which will always == 1.

if false

    %R = R';
    imagesc(abs(R), [0 1]); colormap(cm);
    set(gca, 'TickLength', [0 0]);
    set(gca, 'Xtick', 1:size(R,2));
    set(gca, 'XTickLabel', results.nsamples);

    hold on
    clr = 'g'; %3.3*[0.3 0.3 0.3];
    for i = 1 : size(R,1)
        for j = 1 : size(R,2)
            if R(i,j) > 0.999
                l = j - 0.5;
                r = j + 0.5;
                t = i - 0.5;
                b = i + 0.5;

                plot([l l], [b t], '-', 'Color', clr, 'LineWidth', 3);
                plot([r r], [b t], '-', 'Color', clr, 'LineWidth', 3);
                plot([l r], [t t], '-', 'Color', clr, 'LineWidth', 3);
                plot([l r], [b b], '-', 'Color', clr, 'LineWidth', 3);
            end
        end
    end
    axis image
    %colorbar;
    set(gca, 'FontSize', 14);
    rotateXLabels_imagesc(gca, 45);
    set(gca, 'YTick', 1:NMARKERS);
    set(gca, 'YTickLabel', tids(results.somps{end}.S(1:NMARKERS)));

    plotSave('figures/marker_convergence/marker_convergence_vs_db_size.png');
    close
end


if false
    figure
    colorbar;
    axis off
    colormap(cm);
    caxis([0 1]);
    set(gca, 'FontSize', 14);
    plotSave('figures/marker_convergence/marker_convergence_colorbar.png');
    close
end

if true
    figure
    b = bar(fliplr(100*diff(1 - [1 results.somps{end}.punexp(1:NMARKERS)])),1);
    set(gca, 'FontSize', 20);
    set(gca, 'XTick', []);
    set(b, 'FaceColor', [0.5 0.5 0.5]);
    set(b, 'EdgeColor', 'w');
    set(gca, 'TickLength', [ 0 0 ]);
    axis tight;
    axis off
    box off
    plotSave('figures/marker_convergence/percent_var_explained.png');
    close
   
end




end

