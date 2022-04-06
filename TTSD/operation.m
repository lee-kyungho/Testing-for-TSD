
function f = operation(n,m,sample,grid)

% Operator for calculating EDFT
% Author: Kyungho Lee, Oliver Linton, and Yoon-Jae Whang
% n : Time order
% m : Stochastic order
% Sample
% grid: grid points

%---- About Input Sample ----%
% The input sample should have a dimension as N_k * (T+1) * b

T = size(sample,2)-1; % Terminal period T
b = size(sample,3);   % Bootstrap sample size

% operator
% Note that if n = 1 -> (t+1-s)^(n-1) = 1
% Note that if n = 2 -> (t+1-s)^(n-1) = (t+1-s)

linear_operator = @(x,t,s) (sample(:,1:t+1,:)<=x).*((x-sample(:,1:t+1,:)).^(m-1)).*((t+1-s).^(n-1))/(factorial(m-1));         


for t = 0:T
    s = 0:t;
    % let t and s be given
    linear_operator_given_t_s = @(x) linear_operator(x,t,s);
    %disp("lin op")
    %size(linear_operator_given_t_s)
    % apply along grid (apply operator on each grid)
    applying_operator_cells = arrayfun(linear_operator_given_t_s,grid,'UniformOutput',false);
    %disp("app along grid")
    %size(applying_operator_cells{1})
    
    % concatentae each cell to apply mean and sum functions
    applying_operator_array = cat(4,applying_operator_cells{:});
    %disp("concat grid")
    %size(applying_operator_array)
        
    ecdf_t = nanmean(applying_operator_array,1); % to ignore nan, use nanmean
    %size(ecdf_t)

    time_cumulated_ecdf_t = sum(ecdf_t,2);
    %size(time_cumulated_ecdf_t)

    operation_result(:,t+1,:,:) = time_cumulated_ecdf_t; 
    
end
% Output dim: n_grid x T+1
% Output dim: b x T+1 x n_grid

if b == 1
    f = operation_result; % Output dim: n_grid x T+1
else
    f = permute(operation_result,[3,2,1]); % Output dim: n_grid x T+1 x b
end
end
