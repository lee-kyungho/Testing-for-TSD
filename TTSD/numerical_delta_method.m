
function NDM_dist = numerical_delta_method(D,b_D_recentered,epsilon,r_N,p,type,way)
% Author: Kyungho Lee, Oliver Linton, and Yoon-Jae Whang
% D : difference of EDFTs
% b_D_recentered : recentered differences (of bootstrap samples)
% epsilon : length parameter for NDM
% r_N : sqrt(N)
% p : 1 or 2 (determining the order of L_p statistic)
% type : 'max' or 'sum'
% way : 1 or 2 (the way to approximate the second-order NDM. This parameter does not affect the result when NDM is first-order approximation i.e., if p=1)
% input data dimension: (# of GRID) x (Length of Time Horizon + 1) x (# of Bootstrap Sample)
% If there is no bootstrap sample (i.e. original sample) then the last
% dimension is 1.

% Specify the functional phi
if type == 'max'
phi = @(theta) trapz(Lambda(theta,p,'max'),1);
elseif type == 'sum'
phi = @(theta) trapz(Lambda(theta,p,'sum'),1);
end

phi_D = phi(D); % 1 x 1 

% Output:
% 1 x 1 x (# of Bootstrap sample)
if p == 1

NDM_dist = (phi(D + epsilon*r_N*(b_D_recentered))-phi_D)/epsilon; % Use recentered b_D

elseif p == 2 % higher-order

% Output:
% 1 x 1 x (# of Bootstrap sample)

    if way == 1
        NDM_dist = (phi(D + epsilon*r_N*b_D_recentered) - phi_D) / (epsilon^2);
    elseif way == 2
        NDM_dist = (phi(D + 2*epsilon*r_N*b_D_recentered) - 2*phi(D + epsilon*r_N*b_D_recentered) + phi_D) / 2*(epsilon^2);
    end
end
