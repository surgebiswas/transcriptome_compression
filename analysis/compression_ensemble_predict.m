function yhat = compression_ensemble_predict( xq, model )


% Re-label indices using the idx_all as a basis.
% 0 = anchor indices.
% 1 ... n = context specific compression indices
idx_rl = -ones(1, length(model.idx_all));
idx_rl(1:length(model.idx_anchor)) = 0;
itr = length(model.idx_anchor) + 1;
for i = 1 : length(model.idx_cs)
    idx_rl(itr:itr+length(model.idx_cs{i})-1) = i;
    itr = itr+length(model.idx_cs{i});
end



%%% Anchor reconstruction
b_anchor = model.b_anchor;
yhat_a = xq(:,idx_rl==0)*b_anchor; % Assume anchors are arranged first.
r = xq - yhat_a(:,model.idx_all);

%%% Full reconstruction
a = model.alpha*model.alpha_scale; % Kernel (defined as precision)

% Create weight matrix
w = zeros(size(xq,1),length(model.b_cs));
for i = 1 : length(model.b_cs)
    cs_centered = -bsxfun(@minus, r, model.means(i,:));
    for k = 1 : size(cs_centered,1)
        w(k,i) = -cs_centered(k,:)*a*(cs_centered(k,:)'); % squared Mahalanobis distance
    end
end

% now subtract off the max 
for i = 1 : size(w,1)
    w(i,:) = smax(w(i,:));
end


% Construct predictions.
yhat = yhat_a;
for j = 1 : length(model.b_cs) % over compressions
    yhat = yhat + bsxfun(@times, r(:,idx_rl==j)*model.b_cs{j}, w(:,j));
end



end

