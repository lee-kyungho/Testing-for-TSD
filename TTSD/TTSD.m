% MATLAB codes for implementing 'Testing for Time Stochastic Dominance'
% Author: Kyungho Lee, Oliver Linton, and Yoon-Jae Whang

T = 4;
ngrid = 100;
btsp = 200; % Number of bootstrapiing
p1=1;  %  L1 staistics
p2=2;  %  L2 staistics

% Significance Level
alpha = 0.05;

% eta
eta = 10^(-6);

N = 100;

% Generating Data
% Mean-Shifting Uniform Distributions (as in Section 7)
X = zeros(N,T+1);
X0 = unifrnd(0,1,N,1);
X(:,1,:) = X0;

for i = 1:T
    X(:,i+1,:) = 0.5*X(:,i,:) +  unifrnd(0,1,N,1)  - 1/2;
end

Y = zeros(N,T+1);
Y0 = unifrnd(0,1,N,1);
Y(:,1,:) = Y0;

for i = 1:T
    Y(:,i+1,:) = 0.5*Y(:,i,:) + unifrnd(0,1,N,1) - 1/2;
end

% ------ Testing Starts ------ %

r_N = sqrt(size(X,1));

% Samples
sample1 = squeeze(X(:,:));
sample2 = squeeze(Y(:,:));

% Setting a grid
grid = linspace(min(min(sample1,sample2),[],'all'),max(max(sample1,sample2),[],'all'),ngrid)';
% Function for linear operator

% Bootstrapping 'Path-Wise'
b1sample1 = path_wise_bootstrap(sample1,btsp);
b2sample2 = path_wise_bootstrap(sample2,btsp);

% Calculate Test Statistics

% (Time, SD) order = (1,1)
op1_11 = operation(1,1,sample1,grid); % Output dim ngrid * (T+1) * btsp
op2_11 = operation(1,1,sample2,grid); % Output dim ngrid * (T+1) * btsp
D_11 = op1_11 - op2_11;

% (Time, SD) order = (2,2)
op1_22 = operation(2,2,sample1,grid);
op2_22 = operation(2,2,sample2,grid);
op1_12_T = operation_T(1,2,sample1,grid);
op2_12_T = operation_T(1,2,sample2,grid);

D_22 = op1_22 - op2_22;
D_12_T = op1_12_T - op2_12_T;

% p = 1 (p1)
Lamb_11_max_p1 = Lambda(D_11,p1,'max'); % Output dim: ngrid * 1
% For n=2, m=2, we need to make collection of D_22 and D_12_T

D_22_collection = cat(2,D_22,D_12_T); % ngrid * J

Lamb_22_max_p1 = Lambda(D_22_collection,p1,'max'); 

%---- Integration for T_{N} ----------%

T_11_max_p1 = r_N^p1 * trapz(Lamb_11_max_p1);
T_22_max_p1 = r_N^p1 * trapz(Lamb_22_max_p1);

%---- Bootstrap Sample Statistics -----%

b_op1_11 = operation(1,1,b1sample1,grid);
b_op2_11 = operation(1,1,b2sample2,grid);
b_D_11 = b_op1_11 - b_op2_11;

b_op1_22 = operation(2,2,b1sample1,grid);
b_op2_22 = operation(2,2,b2sample2,grid);
b_D_22 = b_op1_22 - b_op2_22;

b_op1_12_T = operation_T(1,2,b1sample1,grid);
b_op2_12_T = operation_T(1,2,b2sample2,grid);
b_D_12_T = b_op1_12_T - b_op2_12_T;

%----------- Recentering ------------%

b_D_11_recentered = b_D_11 - D_11; 
b_D_22_recentered = b_D_22 - D_22; 
b_D_12_T_recentered = b_D_12_T - D_12_T; 

%---------- Calculating Lambda ------%

% Bootstrap Sample ------------------%

% p = 1

btsp_Lamb_11_max_LFC_p1 = Lambda(b_D_11_recentered,p1,'max'); % Output dim: ngrid * 1 * btsp
b_D_12_T_recentered = reshape(b_D_12_T_recentered,[size(b_D_12_T_recentered,1),1,size(b_D_12_T_recentered,2)]);
b_D_22_collection_recentered = cat(2,b_D_22_recentered,b_D_12_T_recentered); % ngrid * J * btsp
btsp_Lamb_22_max_LFC_p1 = Lambda(b_D_22_collection_recentered,p1,'max'); % Output dim: ngrid * 1 * btsp

%---- Integration for Bootstrap version of T_{N} %
% Bootstrap Sample
% p = 1

% -------- LFC Approach -------- %
btsp_T_11_max_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_11_max_LFC_p1,1);
btsp_T_22_max_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_22_max_LFC_p1,1);

% -------- Contact-set Approach -------- %
c_N_1 = 0.5*log(log(N));

[T_11_max_contact_1_p1] = ...
            TTSD_contact(D_11,b_D_11_recentered, p1, N, r_N, c_N_1, "max");

[T_22_max_contact_1_p1] = ...
            TTSD_contact(D_22,b_D_22_collection_recentered, p1, N, r_N, c_N_1, "max");

% -------- Numerical Delta Method -------- %
epsilon_1 = r_N^(-1/128);

btsp_phi_dist_11_max_1_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_1,r_N,p1,'max',1);
btsp_phi_dist_22_max_1_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_1,r_N,p1,'max',1);

% -------- Calculating Critical Value -------- %

critical_value_11_btsp_max_LFC_p1       = quantile(btsp_T_11_max_LFC_p1,1-alpha);
critical_value_22_btsp_max_LFC_p1       = quantile(btsp_T_22_max_LFC_p1,1-alpha);

critical_value_11_max_contact_1_p1      = max(quantile(T_11_max_contact_1_p1,1-alpha),eta);
critical_value_22_max_contact_1_p1      = max(quantile(T_22_max_contact_1_p1,1-alpha),eta);

critical_value_11_btsp_max_NDM_1_p1     = quantile(btsp_phi_dist_11_max_1_p1,1-alpha);
critical_value_22_btsp_max_NDM_1_p1     = quantile(btsp_phi_dist_22_max_1_p1,1-alpha);

% -------- Deciding Rejection of H0 -------- %

rejection_11_btsp_max_LFC_p1     = T_11_max_p1 > critical_value_11_btsp_max_LFC_p1;
rejection_22_btsp_max_LFC_p1     = T_22_max_p1 > critical_value_22_btsp_max_LFC_p1;

rejection_11_max_contact_1_p1    = T_11_max_p1 > critical_value_11_max_contact_1_p1;
rejection_22_max_contact_1_p1    = T_22_max_p1 > critical_value_22_max_contact_1_p1;

rejection_11_btsp_max_NDM_1_p1   = T_11_max_p1 > critical_value_11_btsp_max_NDM_1_p1;
rejection_22_btsp_max_NDM_1_p1   = T_22_max_p1 > critical_value_22_btsp_max_NDM_1_p1;

% -------- P - value -------- %

% LFC
p_value_11_max_LFC = mean(T_11_max_p1 <= btsp_T_11_max_LFC_p1);
p_value_22_max_LFC = mean(T_22_max_p1 <= btsp_T_22_max_LFC_p1);

% Contact
p_value_11_max_contact = mean(T_11_max_p1 <= T_11_max_contact_1_p1);
p_value_22_max_contact = mean(T_22_max_p1 <= T_22_max_contact_1_p1);

% NDM
p_value_11_max_NDM = mean(T_11_max_p1 <= btsp_phi_dist_11_max_1_p1);
p_value_22_max_NDM = mean(T_22_max_p1 <= btsp_phi_dist_22_max_1_p1);

% -------- Print Results -------- %

fprintf("Time stochastic dominance testing results (%d, %d) order \n", 1, 1)
disp("Procedure :    LFC    Contact    NDM")
fprintf("P-value   :    %.3f    %.3f    %.3f \n",p_value_11_max_LFC, p_value_11_max_contact, p_value_11_max_NDM)
fprintf("Reject    :    %d    %d    %d \n", rejection_11_btsp_max_LFC_p1, rejection_11_max_contact_1_p1, rejection_11_btsp_max_NDM_1_p1)


fprintf("Time stochastic dominance testing results (%d, %d) order \n", 2, 2)
disp("Procedure :    LFC    Contact    NDM")
fprintf("P-value   :    %.3f    %.3f    %.3f \n",p_value_22_max_LFC, p_value_22_max_contact, p_value_22_max_NDM)
fprintf("Reject    :    %d    %d    %d \n", rejection_22_btsp_max_LFC_p1, rejection_22_max_contact_1_p1, rejection_22_btsp_max_NDM_1_p1)
