
clear
% Author: Kyungho Lee, Oliver Linton, Yoon-Jae Whang
% Date: 2020-11-03
% This code is to apply TSD_Testing to Thaipanel data (Kaboski and Townsend 2012, AEJ)

Thai_panel = readtable("ThaiVillage_Panel.csv"); % TWO SAMPLES

Thai_panel.id = Thai_panel.case_id;
Thai_panel = Thai_panel(Thai_panel.maleh ~= 0.5, :);

% Dropping NaN values
Thai_mat = table2array(Thai_panel(:,['netinc','tc',"year",'madult', "fadult", 'kids', 'maleh', 'farm', 'ageh', 'age2h', 'educh']));
Nan_index = logical(sum(isnan(Thai_mat),2));
Thai_panel = Thai_panel(~Nan_index,:);

% Index for Before/After Program
after_program_mask = (Thai_panel.year >= 2002);
after_program_mask = logical(after_program_mask);
before_program_mask = (Thai_panel.year < 2002);
before_program_mask = logical(before_program_mask);


%-------------------------------% 

% TSD testing Starts

%-------------------------------%

for after = [1,0] % Whether we are going to test 'after' the program

if after == 1
    
    Thai_data = Thai_panel(after_program_mask,:);
    year_dummy = dummyvar(categorical(Thai_data.year));
    Thai_data.y2002 = year_dummy(:,1);
    Thai_data.y2003 = year_dummy(:,2);
    Thai_data.y2004 = year_dummy(:,3);
    Thai_data.y2005 = year_dummy(:,4);
    Thai_data.y2006 = year_dummy(:,5);
    Thai_data.y2007 = year_dummy(:,6);

else
    
    Thai_data = Thai_panel(before_program_mask,:);
    year_dummy = dummyvar(categorical(Thai_data.year));
    Thai_data.y1997 = year_dummy(:,1);
    Thai_data.y1998 = year_dummy(:,2);
    Thai_data.y1999 = year_dummy(:,3);
    Thai_data.y2000 = year_dummy(:,4);
    Thai_data.y2001 = year_dummy(:,5);

end


Results = ["design" 'variable' 'order' 'contact'];

for residual_flag = [1]    
    
for variable = ["netinc",  "tc"] % "tc" , 

small_mask = Thai_data.small == 1;
big_mask = Thai_data.small == 0;

Thai_small = Thai_data(small_mask,:);
Thai_big = Thai_data(big_mask,:);

[y_small_Thai, small_unique_hh] = construct_test_data(Thai_small,variable);
small_id = [[1:size(small_unique_hh,1)]', small_unique_hh];
[y_big_Thai, big_unique_hh] = construct_test_data(Thai_big,variable);
big_id = [[1:size(big_unique_hh,1)]', big_unique_hh];

if residual_flag == 1
    
% Linear Regression
if after == 1

Thai_temp = Thai_data(:,[variable, 'madult', "fadult", 'kids', 'maleh', 'farm', 'ageh', 'age2h', 'educh',...
            "y2002","y2003","y2004","y2005","y2006","y2007"]);
Thai_array = table2array(Thai_temp);

y = Thai_array(:,1);
X = Thai_array(:,2:end);

beta_hat = (X'*X)\X'*y;

Thai_temp_small = Thai_small(:,[variable, 'madult', "fadult", 'kids', 'maleh', 'farm', 'ageh', 'age2h', 'educh',...
            "y2002","y2003","y2004","y2005","y2006","y2007"]);
Thai_temp_big   = Thai_big(:,[variable, 'madult', "fadult", 'kids', 'maleh', 'farm', 'ageh', 'age2h', 'educh',...
            "y2002","y2003","y2004","y2005","y2006","y2007"]);

Thai_array_small = table2array(Thai_temp_small);
Thai_array_big = table2array(Thai_temp_big);

variable_y = variable;
variable = "res_" + variable_y;

Thai_small{:,variable} =  Thai_array_small(:,1) - Thai_array_small(:,2:end)*beta_hat;
Thai_big{:,variable} =  Thai_array_big(:,1) - Thai_array_big(:,2:end)*beta_hat;

else 

Thai_temp = Thai_data(:,[variable, 'madult', "fadult", 'kids', 'maleh', 'farm', 'ageh', 'age2h', 'educh',...
            "y1997","y1998","y1999","y2000","y2001"]);
Thai_array = table2array(Thai_temp);

y = Thai_array(:,1);
X = Thai_array(:,2:end);

beta_hat = (X'*X)\X'*y;    

Thai_temp_small = Thai_small(:,[variable, 'madult', "fadult", 'kids', 'maleh', 'farm', 'ageh', 'age2h', 'educh',...
            "y1997","y1998","y1999","y2000","y2001"]);
Thai_temp_big   = Thai_big(:,[variable, 'madult', "fadult", 'kids', 'maleh', 'farm', 'ageh', 'age2h', 'educh',...
            "y1997","y1998","y1999","y2000","y2001"]);

Thai_array_small = table2array(Thai_temp_small);
Thai_array_big = table2array(Thai_temp_big);

variable_y = variable;
variable = "res_" + variable_y;

Thai_small{:,variable} =  Thai_array_small(:,1) - Thai_array_small(:,2:end)*beta_hat;
Thai_big{:,variable}   =  Thai_array_big(:,1) - Thai_array_big(:,2:end)*beta_hat;

end

% Reshaping data set for path-wise bootstrap later

[y_small_Thai, ~] = construct_test_data(Thai_small,variable_y);
[y_big_Thai, ~] = construct_test_data(Thai_big,variable_y);

[res_small_Thai, ~]  = construct_test_data(Thai_small,variable);
[res_big_Thai, ~] = construct_test_data(Thai_big,variable);

[year_small_Thai, ~]  = construct_test_data(Thai_small,"year");
[year_big_Thai, ~] = construct_test_data(Thai_big,"year");

[madult_small_Thai, ~]  = construct_test_data(Thai_small,"madult");
[madult_big_Thai, ~] = construct_test_data(Thai_big,"madult");

[fadult_small_Thai, ~]  = construct_test_data(Thai_small,"fadult");
[fadult_big_Thai, ~] = construct_test_data(Thai_big,"fadult");

[kids_small_Thai, ~]  = construct_test_data(Thai_small,"kids");
[kids_big_Thai, ~] = construct_test_data(Thai_big,"kids");

[maleh_small_Thai, ~]  = construct_test_data(Thai_small,"maleh");
[maleh_big_Thai, ~] = construct_test_data(Thai_big,"maleh");

[farm_small_Thai, ~]  = construct_test_data(Thai_small,"farm");
[farm_big_Thai, ~] = construct_test_data(Thai_big,"farm");

[ageh_small_Thai, ~]  = construct_test_data(Thai_small,"ageh");
[ageh_big_Thai, ~] = construct_test_data(Thai_big,"ageh");

[age2h_small_Thai, ~]  = construct_test_data(Thai_small,"age2h");
[age2h_big_Thai, ~] = construct_test_data(Thai_big,"age2h");

[educh_small_Thai, ~]  = construct_test_data(Thai_small,"educh");
[educh_big_Thai, ~] = construct_test_data(Thai_big,"educh");


% Container of covariates

covariate_small_cg = cat(3,year_small_Thai);
covariate_small = cat(3,  madult_small_Thai, fadult_small_Thai, ...
                         kids_small_Thai, maleh_small_Thai, farm_small_Thai, ageh_small_Thai, age2h_small_Thai,educh_small_Thai);

covariate_big_cg = cat(3,year_big_Thai);
covariate_big = cat(3,  madult_big_Thai, fadult_big_Thai, ...
                         kids_big_Thai, maleh_big_Thai, farm_big_Thai, ageh_big_Thai, age2h_big_Thai, educh_big_Thai);

end

% Tuning Parameters
p = 1; % L1-Type stat

% Bootstrapping
btsp = 1;

for design = [1,2]

disp(design)    

if design == 1
    
    % small >tsd big
    
if residual_flag == 1
    sample1 = squeeze(res_small_Thai(:,:));
    sample2 = squeeze(res_big_Thai(:,:));
    sample1_id = small_id;
    sample2_id = big_id;
    
    sample1_y = y_small_Thai;
    sample2_y = y_big_Thai;

    sample1_covariate_cg = covariate_small_cg; 
    sample1_covariate = covariate_small;
    sample2_covariate_cg = covariate_big_cg;
    sample2_covariate = covariate_big;

else
    sample1 = squeeze(y_small_Thai(:,:));
    sample2 = squeeze(y_big_Thai(:,:));

end

elseif design == 2
    
    % big >tsd small
if residual_flag == 1
    sample1 = squeeze(res_big_Thai(:,:));
    sample2 = squeeze(res_small_Thai(:,:));
    sample1_id = big_id;
    sample2_id = small_id;
    
    sample1_y = y_big_Thai;
    sample2_y = y_small_Thai;

    sample1_covariate_cg = covariate_big_cg;
    sample1_covariate = covariate_big;
    sample2_covariate_cg = covariate_small_cg;
    sample2_covariate = covariate_small;

else
    sample1 = squeeze(y_big_Thai(:,:));
    sample2 = squeeze(y_small_Thai(:,:));
end

end

% Number of grid points
ngrid = 100;

% Grid Points
grid = linspace(min(min(min(sample1),min(sample2)),[],'all'),max(max(max(sample1),max(sample2)),[],'all'),ngrid)';

% Estimation Starts

T = size(sample1,2)-1;
N1 = size(sample1,1);
N2 = size(sample2,1);

r_N = sqrt(N1*N2/(N1+N2));

T_11_sum_contact_2_list = [];
T_12_sum_contact_2_list = [];
T_21_sum_contact_2_list = [];
T_22_sum_contact_2_list = [];

% (Time, SD) order = (1,1)
op1_11 = operation(1,1,sample1,grid); % Output dim ngrid * (T+1) * b
op2_11 = operation(1,1,sample2,grid); % Output dim ngrid * (T+1) * b
D_11 = op1_11 - op2_11;

% (Time, SD) order = (1,2)
op1_12 = operation(1,2,sample1,grid);
op2_12 = operation(1,2,sample2,grid);
D_12 = op1_12 - op2_12;

% (Time, SD) order = (2,1)
op1_21 = operation(2,1,sample1,grid);
op2_21 = operation(2,1,sample2,grid);
D_21 = op1_21 - op2_21;

% (Time, SD) order = (2,2)
op1_22 = operation(2,2,sample1,grid);
op2_22 = operation(2,2,sample2,grid);
op1_12_T = operation_T(1,2,sample1,grid);
op2_12_T = operation_T(1,2,sample2,grid);

D_22 = op1_22 - op2_22;
D_12_T = op1_12_T - op2_12_T;

D_22_collection = cat(2,D_22,D_12_T); % ngrid * J

Lamb_11_sum = Lambda(D_11,p,'sum'); % Output dim: ngrid * 1
Lamb_12_sum = Lambda(D_12,p,'sum'); % Output dim: ngrid * 1
Lamb_21_sum = Lambda(D_21,p,'sum'); % Output dim: ngrid * 1

% For n=2, m=2, we need to make collection of D_22 and D_12_T
Lamb_22_sum = Lambda(D_22_collection,p,'sum'); % Output dim: ngrid * 1

%---- Integration for T_{N} ----------%

T_11_sum = r_N^p * trapz(Lamb_11_sum);
T_12_sum = r_N^p * trapz(Lamb_12_sum);
T_21_sum = r_N^p * trapz(Lamb_21_sum);
T_22_sum = r_N^p * trapz(Lamb_22_sum);
    

parfor b = 1:200

% Bootstrapping 'Path-Wise' 
index_for_path_wise_construnction1 = 0:N1:N1*T;
bootstrap_index1 = randi([1,N1],N1,1,btsp);
pw_btsp_index1 = bootstrap_index1 + repmat(index_for_path_wise_construnction1,1,1,btsp);

index_for_path_wise_construnction2 = 0:N2:N2*T;
bootstrap_index2 = randi([1,N2],N2,1,btsp);
pw_btsp_index2 = bootstrap_index2 + repmat(index_for_path_wise_construnction2,1,1,btsp);

% Check whether we have to conduct regression

if residual_flag == 1

% sample 1
y_bsample1 = sample1_y(pw_btsp_index1);
covr_bsample1 = sample1_covariate(bootstrap_index1,:,:);
covr_cg_bsample1 = sample1_covariate_cg(bootstrap_index1,:,:);

% Unstacking for regression
y_b1 = reshape(y_bsample1,[size(y_bsample1,1)*size(y_bsample1,2),1]);
covr_b1 = reshape(covr_bsample1,[size(covr_bsample1,1)*size(covr_bsample1,2),size(covr_bsample1,3)]);
covr_cg_b1 = reshape(covr_cg_bsample1,[size(covr_cg_bsample1,1)*size(covr_cg_bsample1,2),size(covr_cg_bsample1,3)]);

% sample 2
y_bsample2 = sample2_y(pw_btsp_index2);
covr_bsample2 = sample2_covariate(bootstrap_index2,:,:);
covr_cg_bsample2 = sample2_covariate_cg(bootstrap_index2,:,:);

% Unstacking for regression
y_b2 = reshape(y_bsample2,[size(y_bsample2,1)*size(y_bsample2,2),1]);
covr_b2 = reshape(covr_bsample2,[size(covr_bsample2,1)*size(covr_bsample2,2),size(covr_bsample2,3)]);
covr_cg_b2 = reshape(covr_cg_bsample2,[size(covr_cg_bsample2,1)*size(covr_cg_bsample2,2),size(covr_cg_bsample2,3)]);
 
% OLS Regression Starts
y_b = [y_b1; y_b2];
covr_b = [covr_b1; covr_b2];
covr_cg_b = [covr_cg_b1; covr_cg_b2];

year_dummy_b = dummyvar(categorical(covr_cg_b));
Thai_array_b = [y_b, covr_b, year_dummy_b];

y_b = Thai_array_b(:,1);
X_b = Thai_array_b(:,2:end);

beta_hat_b = pinv(X_b'*X_b)*X_b'*y_b;

X_b1 = X_b(1:size(y_b1,1),:);
X_b2 = X_b(size(y_b1,1)+1:end,:);

b1sample1      = y_b1 - X_b1*beta_hat_b;
b2sample2      = y_b2 - X_b2*beta_hat_b;

b1sample1 = b1sample1(pw_btsp_index1);
b2sample2 = b2sample2(pw_btsp_index2);

else

b1sample1 = sample1(pw_btsp_index1);
b2sample2 = sample2(pw_btsp_index2);

end

%---- Bootstrap Sample Statistics -----%
%---- The Least Favorable Case   -----%

b_op1_11 = operation(1,1,b1sample1,grid);
b_op2_11 = operation(1,1,b2sample2,grid);
b_D_11 = b_op1_11 - b_op2_11;

b_op1_12 = operation(1,2,b1sample1,grid);
b_op2_12 = operation(1,2,b2sample2,grid);
b_D_12 = b_op1_12 - b_op2_12;

b_op1_21 = operation(2,1,b1sample1,grid);
b_op2_21 = operation(2,1,b2sample2,grid);
b_D_21 = b_op1_21 - b_op2_21;

b_op1_22 = operation(2,2,b1sample1,grid);
b_op2_22 = operation(2,2,b2sample2,grid);
b_D_22 = b_op1_22 - b_op2_22;

b_op1_12_T = operation_T(1,2,b1sample1,grid);
b_op2_12_T = operation_T(1,2,b2sample2,grid);
b_D_12_T = b_op1_12_T - b_op2_12_T;

%----------- Recentering ------------%

b_D_11_recentered = b_D_11 - D_11; 
b_D_12_recentered = b_D_12 - D_12; 
b_D_21_recentered = b_D_21 - D_21; 
b_D_22_recentered = b_D_22 - D_22; 
b_D_12_T_recentered = b_D_12_T - D_12_T; 

%---------- Calculating Lambda ------%
Lamb_11_sum_LFC = Lambda(b_D_11_recentered,p,'sum'); % Output dim: ngrid * 1 * b
Lamb_12_sum_LFC = Lambda(b_D_12_recentered,p,'sum'); % Output dim: ngrid * 1 * b
Lamb_21_sum_LFC = Lambda(b_D_21_recentered,p,'sum'); % Output dim: ngrid * 1 * b

% For n = 2, m = 2, we need to make collection of D_22 and D_12_T
% Reshaping for concatenation
b_D_12_T_recentered = reshape(b_D_12_T_recentered,[size(b_D_12_T_recentered,1),1,size(b_D_12_T_recentered,2)]);
b_D_22_collection_recentered = cat(2,b_D_22_recentered,b_D_12_T_recentered); % ngrid * J * b

Lamb_22_sum_LFC = Lambda(b_D_22_collection_recentered,p,'sum'); % Output dim: ngrid * 1 * b

%---- Bootsrap Sample Statistics -----%
%---- Contact Set Approach -----------%

N = (N1 + N2) / 2;

% Setting Tuning parameter 
% c = 0.5
c_N_2 = 0.5*log(log(N));

[T_11_sum_contact_2,T_12_sum_contact_2, ...
 T_21_sum_contact_2,T_22_sum_contact_2] = ...
            Thai_contact_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p, N, r_N, c_N_2);

T_11_sum_contact_2_list = [T_11_sum_contact_2_list T_11_sum_contact_2];
T_12_sum_contact_2_list = [T_12_sum_contact_2_list T_12_sum_contact_2];
T_21_sum_contact_2_list = [T_21_sum_contact_2_list T_21_sum_contact_2];
T_22_sum_contact_2_list = [T_22_sum_contact_2_list T_22_sum_contact_2];

end

% P-value 
p_value_11_sum_contact_2 = mean(T_11_sum <= T_11_sum_contact_2_list);
p_value_12_sum_contact_2 = mean(T_12_sum <= T_12_sum_contact_2_list);
p_value_21_sum_contact_2 = mean(T_21_sum <= T_21_sum_contact_2_list);
p_value_22_sum_contact_2 = mean(T_22_sum <= T_22_sum_contact_2_list);


Results = [Results; [design variable "(1,1)" p_value_11_sum_contact_2;
                     design variable "(1,2)" p_value_12_sum_contact_2;
                     design variable "(2,1)" p_value_21_sum_contact_2;
                     design variable "(2,2)" p_value_22_sum_contact_2;
                    ]]

end
end
end
writematrix(Results,"P_Values_time_after" + after + "Thai_Micro_20211024.csv")

end
