function qt = NCBI_SRA_Mmusculus_build_and_analyze_query_table( qtfile )

    if false
        qt = read_ncbi_sra_query_table(qtfile);
        
        % There is an abnormal date for sample SRR826224.
        % It was supposedly released on Jan 1, 2020, but 
        % all other samples in its submission were released
        % on June 05, 2013.
        qt({'SRR826224'},:).release_date = {'Jun 05, 2013'};
        qt({'SRR826224'},:).release_date_num = datenum('Jun 05, 2013');
        
        % Link run metadata
        if true
            md = dataset('XLSfile', 'NCBI_SRA_Mmusculus_metadata.xlsx', 'Sheet', 1, 'ReadObsNames', true, 'ReadVarNames', true);

            % Link metadata with query table.
            omd = get(md, 'ObsNames');
            oqt = get(qt, 'ObsNames');
            qt.class = cell(size(qt,1),1);
            for j = 1 : length(qt.class);
                qind = find(strcmpi(oqt{j}, omd));
                if isempty(qind)
                    qt.class{j} = 'unannotated';
                else
                    assert(length(qind) == 1)
                    qt.class{j} = md.class{qind};
                end
            end
        end
        
        
        save([qtfile, '.mat'], 'qt')
    else
        load([qtfile, '.mat']);
    end
    
    


end

