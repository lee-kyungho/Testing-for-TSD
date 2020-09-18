
function [rn_phi_D,phi_dist] = numerical_delta_method(D,b_D_recentered,epsilon,r_N,p,type)

% Author: Kyungho Lee(SNU Econ, kh.lee@snu.ac.kr)
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
phi_dist = (phi(D + epsilon*r_N*(b_D_recentered))-phi_D)/epsilon; % Use recentered b_D

rn_phi_D = r_N*phi_D; % 1 x 1 
