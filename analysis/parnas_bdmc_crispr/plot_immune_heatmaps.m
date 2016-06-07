function [pris, ci] = plot_immune_heatmaps( xd, y, varargin)
% y = samples x processes.
xd.source_name = regexprep(xd.source_name, 'T\d', '');
treats = {'No LPS', 'LPS 2 hours', 'LPS 4 hours', 'LPS 6 hours'};
pris = setParam(varargin, 'pris', cell(1,4));
ci = setParam(varargin, 'ci', []);

if isempty(ci)
    [~,ci] = hclust(y);
end

figure;
for i = length(treats) : -1 : 1
    m = strcmpi(xd.treatment, treats{i});
    
    xdm = xd(m,:);
    ym = y(m,:);
    
    n = xdm.reg_class == 0;
    [ymn, mu, sig] = standardize(ym(n,:));
    ymp = standardize(ym(~n,:), 'mu', mu, 'std', sig);
    xdmp = xdm(~n,:);
%     if i == 4
%         
%         pri = hclust(ymp);
% 
%         kop = xdmp.source_name(pri);
%     else
%         [~,pri] = ismember(kop, xdmp.source_name);
%     end
    nri = 1:12; %hclust(ymn(1:10,:));
    if isempty(pris{i})
        sv = sign(mean(ymp));
        [~,pri] = sort(ymp*sv', 'ascend');
        pris{i} = pri;
    else
        pri = pris{i};
    end
    
    ytoim = [ymn(nri,ci); ymp(pri,ci)];
    subplot(1,length(treats),5 - i);
    imagesc(ytoim, [-3 2]); colormap(parula)
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'TickLength', [0 0]);
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1]) % stretch its width and height
    
%     
%     yneg = ym(n,:);
%     ypos = ym(~n,:);
%     
%     ytoim = [yneg(1:10,:);ypos
end
    


% m6 = strcmpi(xd.treatment, 'LPS 6 hours');
% 
% %ko2ind = setParam(varargin, 'ko2ind', build_KO_sort_indices(xd(m6,:), y(m6,:)));
% ri = setParam(varargin, 'ri', []);
% ci = setParam(varargin, 'ci', []);
% 
% 
% fh = figure;
% kk = [];
% for i = length(treats) : - 1 : 1
%     
%     
%     neg = xd.reg_class == 0 & strcmpi(xd.treatment, treats{i});
%     pos = xd.reg_class == 1 & strcmpi(xd.treatment, treats{i});
%     
%     neg_mean = mean(y(neg,:));
%     neg_std = std(y(neg,:));
%     
%     yneg = y(neg,:); yneg = yneg(1:12,:);
%     ypos = y(pos,:);
%     xdp = xd(pos,:);
%     
%     
%     % Assign sort indices
% %     pi = zeros(size(xdp,1),1);
% %     for p = 1 : size(xdp,1)
% %         pi(p) = ko2ind(xdp.KO{p});
% %     end
% %     [~,sidx] = sort(pi);
% 
% 
%     
%     
%     if i == 4 %&& isempty(ri) && isempty(ci)   
%         if isempty(ri) && isempty(ci); [ri,ci] = hclust(ypos); end
%         kop = xdp.KO(ri,:);
%         [~,sidx] = ismember(kop, xdp.KO);
%         
%         
%         figure;
%         cmap = [0 0 0; 1 0 0; 0 1 0.4; 0.3 0 1; 0.6 0.6 0.6];
%         cvn = zeros(size(yneg,1),1);
%         cvp = zeros(size(ypos,1),1);
%         
%         xdp2 = xdp(sidx,:);
%         cvp(strcmpi(xdp2.pathway, 'canonical')) = 1;
%         cvp(strcmpi(xdp2.pathway, 'ost')) = 2;
%         cvp(strcmpi(xdp2.pathway, 'paf')) = 3;
%         cvp(strcmpi(xdp2.pathway, 'other')) = 4;
%         imagesc([cvn; cvp]); colormap(cmap);
%         set(gca, 'XTick', []);
%         set(gca, 'YTick', []);
%         set(gca, 'TickLength', [0 0]);
%         axis equal 
%         axis tight
%         daspect([1 4 1]);
%         
%         set(0,'CurrentFigure',fh)
%     end
%     
%     [~,sidx] = ismember(kop, xdp.KO);
%     
%     y_toplot = standardize([yneg; ypos(sidx,:)], 'mu', neg_mean, 'std', neg_std);
%     
%     kk = [kk, xdp.KO(sidx)];
%     
%     subplot(1,length(treats),i)
%     imagesc(y_toplot(:,ci), [-2 2]); colormap(parula)
%     
%     hold on
%     plot([0 size(y_toplot,2)+0.5], sum(size(yneg,1)) + 0.5 + [0 0], '-m', 'LineWidth', 2);
%     set(gca, 'XTick', []);
%     set(gca, 'YTick', []);
%     set(gca, 'TickLength', [0 0]);
%     
%     
%     
% end
% 
% kk
% 
%     function ko2ind = build_KO_sort_indices(x,y)
%         k = x.KO;
%         [uk, ia] = unique(k);
%         
%         % Pull out none
%         ni = find(strcmpi(uk, 'none'));
%         uk(ni) = [];
%         
%         
%         ya = [];
%         for j = 1 : length(uk)
%             mask = strcmpi(k, uk{j});
%             ya = [ya; mean(y(mask,:), 1)];
%         end
%         
%         [rr,cc] = hclust(ya);
%         
%         ko2ind = containers.Map;
%         ko2ind('none') = 1;
%         for j = 1 : length(uk)
%             ko2ind(uk{j}) = rr(j)+1;
%         end
%     end


end

