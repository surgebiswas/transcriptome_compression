function [ somp ] = expression_optimize( sY, somp, mean_expression, target_expression, varargin )
% somp = struct output from marker OMP. 
% mean_expression = mean expression of all genes (vector)

corrthresh = setParam(varargin, 'corrthresh', 0.7); 

assert(~isempty(somp.crosscorr), 'Marker OMP needs to be run with ''crosscorr'' set to ''true''');

sY_S = sY(:,somp.S);

R = corr(sY_S, sY);

S = [];
for i = 1 : size(R,1)
    gidx = find(R(i,:) > corrthresh);
    mexp = mean_expression(gidx);
    
    [~,bestidx] = min( abs(mexp - target_expression) );
    S = [S, gidx(bestidx)];
   
end

somp.S_expopt = S;
somp.punexp_expopt = tv(sY - sY(:,somp.S_expopt)*pinv(sY(:,somp.S_expopt))*sY)/tv(sY); 


    function v = tv(x)
        v = sum(var(x));
    end


end

