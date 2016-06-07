function [ ya ] = parnas_batch_correct( xeff, y, crm )

b = estimate_effects(xeff, y);
ya = y - xeff(:,crm)*b(crm,:);


end

