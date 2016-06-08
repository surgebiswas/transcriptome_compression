function [ r, rsq, me, n ] = generate_trprogram_perf3_report( ta, pa, setnames, c, meanexp )


for i = 1 : size(setnames,1)
    rsq(i) = corr(ta(:,i), pa(:,i));
    
    me(i) = mean(meanexp(c(:,i) ~= 0));
    
    n(i) = sum(c(:,i) ~= 0);
end

r = [setnames(:,1), num2cell(rsq)', num2cell(me)', num2cell(n)'];
rsq = rsq';
me = me';
n = n';

disp(corr([rsq, me, log(n)], 'type', 'spearman'));



end

