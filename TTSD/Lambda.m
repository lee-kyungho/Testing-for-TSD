
function L = Lambda(D_collectioon,p,type)

% D_collection should be n_grid * J dimension * b
% About dimension: 
% if (n,m) = (1,1) , (1,2) , (2,1), J = T+1
% if (n,m) = (2,2), J = (T+1) + (n-1) = n + T

max_D_collection = max(0, D_collectioon).^p;

% max_D_collection should be n_grid * J dimension (the same)
if type == 'max'
    L = max(max_D_collection,[],2); % n_grid * 1
elseif type == 'sum'
    L = sum(max_D_collection,2);
end
end
