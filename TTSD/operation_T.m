
function f = operation_T(n,m,sample,grid)

% Operator for the terminal period
% Author: Kyungho Lee, Oliver Linton, and Yoon-Jae Whang
% n : Time order
% m : Stochastic order
% Sample
% grid: grid points

T = size(sample,2)-1; % Terminal period T
b = size(sample,3);

% operator
% Note that if n = 1 -> (t+1-s)^(n-1) = 1
% Note that if n = 2 -> (t+1-s)^(n-1) = (t+1-s)

linear_operator = @(x,t,s) (sample(:,1:t+1,:)<=x).*((x-sample(:,1:t+1,:)).^(m-1)).*((t+1-s).^(n-1))/(factorial(m-1)); 

s = 0:T;
% let t and s be given
linear_operator_given_t_s = @(x) linear_operator(x,T,s);
% apply along grid (apply operator on each grid)
applying_operator_cells = arrayfun(linear_operator_given_t_s,grid,'UniformOutput',false);
% concatentae each cell to apply mean and sum functions
applying_operator_array = cat(4,applying_operator_cells{:});
ecdf_t = nanmean(applying_operator_array,1); % to ignore nan, use nanmean
time_cumulated_ecdf_t = sum(ecdf_t,2);
operation_result = time_cumulated_ecdf_t; % dimension: 1 x 1 x b x n_grid

operation_result = squeeze(operation_result); % dimension: b x n_grid

if b == 1
    f = operation_result;
else
    f = permute(operation_result,[2,1]); % dimension: n_grid x b
end

end
