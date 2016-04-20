function stats = get_geo_metadata_parnas( geoid )

r = urlread(['http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', geoid]);

s = regexp(r, 'knockout: (.*?)<br>', 'tokens');
stats.ko = set_finding(s);

s = regexp(r, 'tissue/cell line: (.*?)<br>', 'tokens');
stats.tcl = set_finding(s);

s = regexp(r, 'condition: (.*?)<br>', 'tokens');
stats.con = set_finding(s);

s = regexp(r, 'treatment: (.*?)<br>', 'tokens');
stats.treat = set_finding(s);


s = regexp(r, '>Source name(.+?)</tr>', 'tokens');
stats.source_name = 'NA';
if ~isempty(s)
    q = regexp(s{1}{1}, '">(.+?)<b', 'tokens');
    if ~isempty(q)
        stats.source_name = q{1}{1};
    end
end
        


    function ff = set_finding(s)
        if ~isempty(s)
            ff = s{1}{1};
        else
            ff = 'NA';
        end
    end






end

