                       
function [T_contact] = ...
          TTSD_contact(D,b_D_recentered, ...
          p, N, r_N, c_N, stat_type)

% D
% b_D_recentered
% p
% r_N
% c_N: Constant for contact-set band
% stat_type: "max" or "sum"
      
 epsilon = 10^(-6);

% Calculating contact-set width
% s_n equals to r_n{v_n - v*_n}/sigma (sigma is 1 in here)

s_n = r_N.*b_D_recentered;
S_n = max(max(max(s_n)),epsilon*sqrt(log(N)));
alpha_n = 0.1/log(N);
c_n = c_N*quantile(S_n,1-alpha_n);

T_contact = 0;

N_J = 1:size(D,2); % N_J_11 = N_J_12 = N_J_21 
N_subset_colletion = PowerSet(N_J);


for i = 1:size(N_subset_colletion,2)

    A = N_subset_colletion{1,i};
    A_complement = setdiff(N_J,A);
    % B_ Contains index for numbers in grid with contact set for given A and
    % If there is no grid points in Contact set, then we take small
    % numbers for the critical value
    
    % Censoring Test Stat
    
    % censoring Bootstrap Statistics
    b_D_A_recentered = b_D_recentered;    
    b_D_A_recentered(:,A_complement,:) = 0;
  
    B = contact_set_estimation(r_N,D,A,c_n); % censoring is in the function
    % -------- Bootstrap Lambda
    
    Lamb_A = Lambda(b_D_A_recentered,p,stat_type); % Output dim: ngrid * 1 * b

    %---- Integration for Bootstrap version of T_{N} %
    % In here, we take integration on "Contact" set
    
    T_contact = T_contact + r_N^p * trapz(Lamb_A(B,:,:),1);
    
end 