function [T_11_sum_contact, T_12_sum_contact, ...
          T_21_sum_contact, T_22_sum_contact] = ...
          Thai_contact_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
          b_D_21_recentered, b_D_22_collection_recentered, p, N, r_N, c_N)


       
% Author: Kyungho Lee, Oliver Linton and Yoon-Jae Whang
% In this code, we implement TGI Estimation of (1,1), (1,2), (2,1), (2,2)
% by once.

epsilon = 10^(-6); 

% This code is to normalize tuning parameters in data-driven way.
% s_n equals to r_n{v_n - v*_n}/sigma (sigma is 1 in here)

s_n_11 = r_N.*b_D_11_recentered;
s_n_12 = r_N.*b_D_12_recentered;
s_n_21 = r_N.*b_D_21_recentered;
s_n_22 = r_N.*b_D_22_collection_recentered;

S_n_11 = max(max(max(s_n_11)),epsilon*sqrt(log(N)));
S_n_12 = max(max(max(s_n_12)),epsilon*sqrt(log(N)));
S_n_21 = max(max(max(s_n_21)),epsilon*sqrt(log(N)));
S_n_22 = max(max(max(s_n_22)),epsilon*sqrt(log(N)));

alpha_n = 0.1/log(N);

c_n_11 = c_N*quantile(S_n_11,1-alpha_n);
c_n_12 = c_N*quantile(S_n_12,1-alpha_n);
c_n_21 = c_N*quantile(S_n_21,1-alpha_n);
c_n_22 = c_N*quantile(S_n_22,1-alpha_n);

% % We now ignore data-driven part of tuning parameter selection
% c_n_11 = c_N;
% c_n_12 = c_N;
% c_n_21 = c_N;
% c_n_22 = c_N;
T_11_sum_contact = 0;
T_12_sum_contact = 0;
T_21_sum_contact = 0;
T_22_sum_contact = 0;


N_J_11 = 1:size(D_11,2); % N_J_11 = N_J_12 = N_J_21 
N_subset_colletion_11 = PowerSet(N_J_11);

N_J_22 = 1:size(D_22_collection,2); 
N_subset_colletion_22 = PowerSet(N_J_22);


for i = 1:size(N_subset_colletion_11,2)

    A = N_subset_colletion_11{1,i};
    A_complement = setdiff(N_J_11,A);
    % B_ Contains index for numbers in grid with contact set for given A
    
    % censoring Bootstrap Statistics
    b_D_11_A_recentered = b_D_11_recentered;    
    b_D_11_A_recentered(:,A_complement,:) = 0;
    
    b_D_12_A_recentered = b_D_12_recentered;   
    b_D_12_A_recentered(:,A_complement,:) = 0;
    
    b_D_21_A_recentered = b_D_21_recentered;       
    b_D_21_A_recentered(:,A_complement,:) = 0;

    B_11 = contact_estimation(r_N,D_11,A,c_n_11); % censoring is in the function
    B_12 = contact_estimation(r_N,D_12,A,c_n_12); 
    B_21 = contact_estimation(r_N,D_21,A,c_n_21); 
           
    % -------- Bootstrap Lambda
    
    Lamb_11_sum_A = Lambda(b_D_11_A_recentered,p,'sum'); % Output dim: ngrid * 1 * b
    Lamb_12_sum_A = Lambda(b_D_12_A_recentered,p,'sum'); % Output dim: ngrid * 1 * b
    Lamb_21_sum_A = Lambda(b_D_21_A_recentered,p,'sum'); % Output dim: ngrid * 1 * b


    %---- Integration for Bootstrap version of T_{N} %
    % In here, we take integration on "Contact" set
    
    T_11_sum_contact = T_11_sum_contact + r_N^p * trapz(Lamb_11_sum_A(B_11,:,:),1);
    T_12_sum_contact = T_12_sum_contact + r_N^p * trapz(Lamb_12_sum_A(B_12,:,:),1);
    T_21_sum_contact = T_21_sum_contact + r_N^p * trapz(Lamb_21_sum_A(B_21,:,:),1);    
    
end 

% For (n,m) = (2,2)

for i = 1:size(N_subset_colletion_22,2)

    A = N_subset_colletion_22{1,i};
    A_complement = setdiff(N_J_22,A);

    % B_ Contains index for numbers in grid with contact set for given A and
    % If there is no grid points in Contact set, then the critical value
    % would be a small number 'eta'.
    B_22 = contact_estimation(r_N,D_22_collection,A,c_n_22); % Use a collection for considering the terminal period
    
    b_D_22_collection_A_recentered = b_D_22_collection_recentered;       
    b_D_22_collection_A_recentered(:,A_complement,:) = 0;
        
    Lamb_22_sum_A = Lambda(b_D_22_collection_A_recentered,p,'sum'); % Output dim: ngrid * 1 * b
    
    %---- Integration for Bootstrap version of T_{N} %
    % In here, we take integration on "Contact" set
        
    T_22_sum_contact = T_22_sum_contact + r_N^p * trapz(Lamb_22_sum_A(B_22,:,:),1);
    
end