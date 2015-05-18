function [ sucessfully_processed, unsucessfully_processed ] = check_srafish_output( srafish_outdir )
% check_srafish_output - Checks the output of srafish and returns a list of
% successfully and unsucessfully processed SRA entries. Here,
% "sucessfully proccessed" means that sailfish finished quantifying all
% transcripts. This does not necessarily mean that the processed sample was
% of high quality.

VERBOSE = true;

old = cd(srafish_outdir);
listing = dir;
listing(1:2) = []; 

sucessfully_processed = [];
unsucessfully_processed = [];
for i = 1 : length(listing)
    if isdir(listing(i).name)
        if was_successfully_processed(listing(i).name);
            sucessfully_processed{end+1} = listing(i).name;
            fprintf('%s -> successful! :]\n', listing(i).name);
        else
            unsucessfully_processed{end+1} = listing(i).name;
            fprintf('%s -> unsuccessful. :[\n', listing(i).name);
        end
    end
end

cd(old);

    function s = was_successfully_processed(d)
        o = cd(d);
        
        % Check to make sure the bias corrected count file from sailfish
        % was output.
        if exist('quant_bias_corrected.sf', 'file')
            q = dir('quant_bias_corrected.sf');
            
            % Make sure the file is reasonably large, and therefore likely
            % properly outputted.
            if q.bytes > 500000
                s = true;
            else
                s = false;
            end
        else
            s = false; 
        end
        
        cd(o);
    end


end

