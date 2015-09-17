function [ mask ] = remove_samples_with_non_prot_coding_contamination( Ypc, Ynpc, sids )
    
    % Load metadata 
    md = dataset('XLSfile', 'NCBI_SRA_Athaliana_run_metadata.xlsx', 'Sheet', 1, 'ReadObsNames', true, 'ReadVarNames', true);
    polyA = cell(size(sids));
    omd = get(md, 'ObsNames');
    for i = 1 : length(sids)
        qind = find(strcmpi(omd, sids{i}));
        if isempty(qind)
            polyA{i} = 'NA';
        else
            polyA{i} = md.polyA_selection{qind};
        end
    end
    
    r = sum(Ynpc,1)./( sum(Ypc,1) + sum(Ynpc,1) );
    
    [rs,sidx] = sort(r, 'descend');
    polyAs = polyA(sidx);
    
    NA = strcmpi(polyAs, 'NA');
    yes = strcmpi(polyAs, 'yes');
    no = strcmpi(polyAs, 'no');
    
    figure
    sf = get_standard_figure_font_sizes;
    hold on
    hb = bar(rs); ydata = get(hb, 'YData');
    plot(find(no), ydata(no), 'or');
   % plot(find(NA), ydata(NA), '.', 'Color', [0.6 0.6 0.6]);
    xlabel('Sample', 'FontSize', sf.axis_labels);
    ylabel('Proportion non-protein coding', 'FontSize', sf.axis_labels);
    
    
    % The above plot suggests any sample with > ~58% non-protein coding is a good idea
    % to remove.
    min_nonpolyA_prop = min(rs(no));
    CUTOFF = min_nonpolyA_prop - 0.03; % cutoff = 3% less than the non-polyA sample with min contamination.
    hold on
    plot([0 length(rs)], [CUTOFF CUTOFF], '-r', 'LineWidth', 2);
    axis tight;
    axis square
    box on
    set(gca, 'FontSize', sf.axis_tick_labels);
    plotSave('figures/quality_filter/non_protein_coding_proportions_and_cutoff.png');
    close
    
    mask = r > CUTOFF;
    
    
    
    
%     fHand = figure;
%     aHand = axes('parent', fHand);
%     hold(aHand, 'on')
%     for i = 1 : length(rs)
%         if strcmpi(polyAs{i}, 'NA')
%             bar(i, rs(i), 'parent', aHand, 'facecolor', [0.6 0.6 0.6]);
%         elseif strcmpi(polyAs{i}, 'yes')
%             bar(i, rs(i), 'parent', aHand, 'facecolor', 'k');
%         elseif strcmpi(polyAs{i}, 'no');
%             bar(i, rs(i), 'parent', aHand, 'facecolor', 'r');
%         end
%     end
    

end

