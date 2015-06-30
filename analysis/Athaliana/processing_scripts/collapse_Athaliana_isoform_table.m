function dc = collapse_isoform_table(d)
    d = sortrows(d, 'obsnames');
    o = regexprep(get(d, 'ObsNames'), '\.\d+', '');
    M = double(d);
    
    
    [uid,~,ic] = unique(o);
    
    Mc = zeros(length(uid),size(M,2));
    
    for i = 1 : length(ic)
        Mc(ic(i),:) = Mc(ic(i),:) + M(i,:);
    end
    
    dc = mat2dataset(Mc, 'ObsNames', uid, 'VarNames', get(d,'VarNames'));
    
    
    
end