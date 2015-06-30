function qt = NCBI_SRA_Mmusculus_build_and_analyze_query_table( qtfile )

    if true
        qt = read_ncbi_sra_query_table(qtfile);
        
        % There is an abnormal date for sample SRR826224.
        % It was supposedly released on Jan 1, 2020, but 
        % all other samples in its submission were released
        % on June 05, 2013.
        qt({'SRR826224'},:).release_date = {'Jun 05, 2013'};
        qt({'SRR826224'},:).release_date_num = datenum('Jun 05, 2013');
        
        save([qtfile, '.mat'], 'qt')
    else
        load([qtfile, '.mat']);
    end
    



end

