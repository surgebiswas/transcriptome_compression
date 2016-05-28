function [ci, m] = build_jbm_heatmaps_1( xd, m, varargin )
% m = samples x features
    rep_to_use = setParam(varargin, 'rep_to_use', 1);
    ci = setParam(varargin, 'use_column_indices', []);
    genos = {'coi1-16', 'npr1-1', 'Col-0', '35S:HopBB1', 'eds16-1'};
    treats = setParam(varargin, 'treatments', {'BTH', 'MeJA', 'Mock'});
    cmap = setParam(varargin, 'colormap', parula);
    stdize = setParam(varargin, 'standardize', true);
    caxis = setParam(varargin, 'caxis', [-3 3]);
    
    mask = steq(xd.geno, genos) & xd.rep == rep_to_use;
    
    xdm = xd(mask,:);
    xdm.time = str2double(xdm.time);
    if stdize
        m = standardize(m(mask,:)); %standardize(m(mask,:)); %standardize( standardize(m(mask, :)')' );
    else
        m = m(mask,:);
    end
    
    
    if isempty(ci)
        [~,ci] = hclust(m);
    end
    
    
    times = [0 1 5 8];
    plotidx = 1;
    for i = 1 : length(treats)
        for j = 1 : length(times)
            msub = [];
            for k = 1 : length(genos)
                if times(j) == 0
                    q = xdm.time == times(j) & strcmpi(xdm.geno, genos{k});
                else
                    q = strcmpi(xdm.treatment, treats{i}) & xdm.time == times(j) & strcmpi(xdm.geno, genos{k});
                end
                msub = ([msub; mean(m(q,:),1)]);
                disp(genos{k});
            end
            
            if length(treats) == 1 || (times(j) ~= 0 || (times(j) == 0 && strcmpi(treats{i}, 'MeJA')))
                subplot(length(treats), length(times),plotidx);
                sub_pos = get(gca,'position'); % get subplot axis position
                set(gca,'position',sub_pos.*[1 1 1.25 1.2*length(treats)/3]) % stretch its width and height
                imagesc((msub(:,ci)), caxis);
                colormap(cmap);
                set(gca, 'XTick', []);
                set(gca, 'YTick', []);
                set(gca, 'TickLength', [0 0]);
                %!title(['time = ', num2str(times(j)), '; treat = ', treats{i}]);
                
                if length(treats) > 1 && times(j) == 0 && strcmpi(treats{i}, 'MeJA')
                    pos = get(gca, 'Position');
                    xoffset = -0.1;
                    pos(1) = pos(1) + xoffset;
                    set(gca, 'Position', pos)
                end
            end
            
            plotidx = plotidx + 1;
        end
    end
    
    

end

