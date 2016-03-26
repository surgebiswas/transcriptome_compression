function [ yhat ] = ridgefit_predict( x0, model )

x0 = [ones(size(x0,1),1), x0];
yhat = x0*model.b;

end

