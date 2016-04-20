function [Y, sids, tids] = quality_filter( s, qfparams, organism, varargin )
%

genplot = setParam(varargin, 'genplot', false);
reportFile = [organism, '_quality_filter_summary_data.mat'];
reportDesc = ['Report for quality_filter.m run on ', organism, '. '];
save(reportFile, 'reportDesc');


Y = s.tpm; 
sids = s.ids;
tids = s.transcript_id;
depth = s.depth;
mapped_ratio = s.mapped_ratio;


% if strcmpi(organism, 'Athaliana')
%     RCTHRESH = 4e6;
%     MRTHRESH = 0.75;
%     CORRCUTOFF = 0.45; % 0.7
%     NZCUTOFF = 0.2; % 0.35
%     TPMCUTOFF = 1;
% elseif strcmpi(organism, 'Mmusculus')
%     RCTHRESH = qfparams.RCTHRESH;
%     MRTHRESH = qfparams.MRTHRESH;
%     CORRCUTOFF = qfparams.CORRCUTOFF; % 0.7
%     NZCUTOFF = qfparams.NZCUTOFF; % 0.35
%     TPMCUTOFF = qfparams.TPMCUTOFF;
% end

RCTHRESH = qfparams.RCTHRESH;
MRTHRESH = qfparams.MRTHRESH;
CORRCUTOFF = qfparams.CORRCUTOFF; % 0.7
NZCUTOFF = qfparams.NZCUTOFF; % 0.35
TPMCUTOFF = qfparams.TPMCUTOFF;


% SAMPLE FILTERING
% 1. Depth and Mapped depth filter
[Y,sids] = filter_by_depth_and_mapping_rate(Y, sids, depth, mapped_ratio, genplot);
fprintf('After mapped pct/depth threshold: Num. samples = %0.0f, Num. features = %0.0f\n', length(sids), length(tids));

% 2. Collapse to a gene table. 
if strcmpi(organism, 'Athaliana')
    % 2.a.
    Yd = collapse_Athaliana_isoform_table(mat2dataset(Y, 'ObsNames', tids, 'VarNames', sids));
    
    % 2.b. Keep only nuclear protein coding genes.
    dk = dataset('file','At_nuclear_protein_coding.txt', ...
        'ReadObsNames', false, 'ReadVarNames', false);
    
    o = get(Yd, 'ObsNames');
    ok = dk.Var1;

    [~, ia] = intersect(o, ok);
    npc = setdiff(1:length(o), ia);
    Yd_npc = Yd(npc,:); % Non protein coding features.
    Yd = Yd(ia,:);
    
    Y = double(Yd);
    Y_npc = double(Yd_npc);
    tids = get(Yd, 'ObsNames');
    sids = get(Yd, 'VarNames');
    
    % 2.c. Remove samples with non-protein coding contamination
    rmmask = remove_samples_with_non_prot_coding_contamination(Y, Y_npc, sids);
    Y(:,rmmask) = [];
    sids(rmmask) = [];
    
elseif strcmpi(organism, 'Mmusculus');
    [tids, Y] = collapse_mouse_isoform_table(tids, Y);
    
end
fprintf('After collapsing isoform table: Num. samples = %0.0f, Num. features = %0.0f\n', length(sids), length(tids));

% 4. Expression filtering
[Y,tids] = expression_filter(Y, tids, TPMCUTOFF);
fprintf('After Expression filtering: Num. samples = %0.0f, Num. features = %0.0f\n', length(sids), length(tids));


% 5. Correlation filtering and filtering based on number of non-zero protein
% coding or lncRNA features.
[Y,sids] = filter_by_corr_and_nz(Y, sids, CORRCUTOFF, NZCUTOFF, genplot);
fprintf('After correlation filtering: Num. samples = %0.0f, Num. features = %0.0f\n', length(sids), length(tids));


% 6. Coefficient of variation for all genes to ensure that everything is
% relatively perturbed.
coov = std(Y,0,2)./mean(Y,2);
save(reportFile, 'coov', '-append');






    %%%%% --- ACCESSORY FUNCTIONS --- %%%%%
    function [y,tids] = expression_filter(y,tids,tpmcutoff)
        k = mean(y,2) > tpmcutoff; % average > cutoff
        
        nz = sum(y > tpmcutoff, 2); % number of samples in which a gene > cutoff
        
        % To be kept a gene must have an average expression greater than
        % the cutoff AND be above the cutoff in 5% of samples. 
        k = k & nz > 0.05*size(y,2); 
        
        y = y(k,:);
        tids = tids(k);
    end
    
    function [y, sids] = filter_by_corr_and_nz(y,sids, corr_cutoff, nz_cutoff, genplot)
        lY = log10(y' + 0.1);
        R = corr(lY'); clear lY; % save memory
        q = mean(R,2);

        z = sum(y == 0)/size(y,1);
        
        if genplot
            sf = get_standard_figure_font_sizes;
            plot(z,q, '.k');
            set(gca, 'FontSize', sf.axis_tick_labels);
            xlabel('Proportion of transcripts with 0 TPM', 'FontSize', sf.axis_labels);
            ylabel('Average corr. with all other samples', 'FontSize', sf.axis_labels);
            axis square
            hold on
            plot([nz_cutoff nz_cutoff], [corr_cutoff 1], '-r', 'LineWidth', 2);
            plot([0 nz_cutoff], [corr_cutoff corr_cutoff], '-r', 'LineWidth', 2);
            
            plotSave('figures/quality_filter/corr_with_other_samples_vs_nz_prop.png');
            close;
        end
        reportDesc = [reportDesc, 'z = Proportion of zero expressed genes per sample. q = Average correlation to other samples. '];
        save(reportFile, 'z', 'q', 'reportDesc', '-append');
        
        k = q' > corr_cutoff & z < nz_cutoff;
        y = y(:,k);
        sids = sids(k);
    end
    
    function [y, sids] = filter_by_depth_and_mapping_rate(y, sids, d, mr, genplot)
        if genplot
            sf = get_standard_figure_font_sizes;
            semilogx((d), 100*(mr), '.k');
            set(gca, 'FontSize', sf.axis_tick_labels);
            xlabel('Mapped depth', 'FontSize', sf.axis_labels);
            ylabel('Mapped percentage', 'FontSize', sf.axis_labels);
            axis square
            hold on
            plot([RCTHRESH RCTHRESH], 100*[MRTHRESH 1], '-r', 'LineWidth', 2);
            plot([RCTHRESH 1e10], 100*[MRTHRESH MRTHRESH], '-r', 'LineWidth', 2);
            
            v = axis;
            v(2) = 2e8;
            v(1) = 100;
            axis(v);
            set(gca, 'XTick', round(logspace(2,8,7)));

            plotSave('figures/quality_filter/mapped_pct_vs_depth.png');
            close;
        end
        reportDesc = [reportDesc, 'd = depth, mr = mapping rate. '];
        save(reportFile, 'd', 'mr', 'reportDesc', '-append');
        
        
        torm = d < RCTHRESH & mr < MRTHRESH;
        y(:, torm) = [];
        sids(torm) = [];
    end



end

