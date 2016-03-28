function [ model ] = treebag_reg_train( x, y, params)
% params{1} = number of trees.

model = TreeBagger(params{1}, x, y, 'method', 'regression');


end

