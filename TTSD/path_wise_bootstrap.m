function b_sample = path_wise_bootstrap(sample,b)

% Function for "path-wise" bootstrap

N = size(sample,1); % sample size of sample
T = size(sample,2) - 1;  % terminal time period of sample1 and sample2

% Bootstrapping 'Path-Wise'
% N = size(sample1,1);
index_for_path_wise_construnction = 0:N:N*T;
boostrap_index = randi([1,N],N,1,b);
boostrap_index = boostrap_index + repmat(index_for_path_wise_construnction,1,1,b);
b_sample = sample(boostrap_index);




