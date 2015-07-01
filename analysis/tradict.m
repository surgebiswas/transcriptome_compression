function yhat = tradict( x, model )

yhat = x*model.b + repmat(model.b0, size(x,1), 1);

end

