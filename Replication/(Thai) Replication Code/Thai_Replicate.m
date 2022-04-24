
% Author: Kyungho Lee, Oliver Linton, Yoon-Jae Whang
% Date: 2020-11-03
% This code is to apply TSD_Testing to Thaipanel data (Kaboski and Townsend 2012, AEJ)

Thai_panel = readtable("ThaiVillage_Panel.csv"); % TWO SAMPLES

% Thai_panel = readtable("ThaiVillage_Panel_Q5.csv");

Thai_panel.id = Thai_panel.case_id;
Thai_panel = Thai_panel(Thai_panel.maleh ~= 0.5, :);

after_program_mask = (Thai_panel.year >= 2002);
after_program_mask = logical(after_program_mask);
before_program_mask = (Thai_panel.year < 2002);
before_program_mask = logical(before_program_mask);

%-------------------------------% 

% TSD testing Starts

%-------------------------------%

for after = [0,1] % Whether we are going to test 'after' the program

if after == 1
    
    Thai_data = Thai_panel(after_program_mask,:);
else
    
    Thai_data = Thai_panel(before_program_mask,:);

end

Thai_data.year = categorical(Thai_data.year);
Results = ["design" 'variable' 'order' 'contact'];

for residual_flag = [0,1]    
    
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

fitted = fitlm(Thai_data, variable + '~ year + madult + fadult + kids + maleh'...
                                + ' + farm + ageh + age2h + educh ' ...                           
                                ,'Intercept',false); % OLS regression 
% Extracting Residuals
                          
variable_y = variable;
variable = "res_" + variable_y;

Thai_data{:,variable} = fitted.Residuals.Raw;

% small_mask = Thai_data.small == 1;
% big_mask = Thai_data.small == 0;

Thai_small = Thai_data(small_mask,:);
Thai_big = Thai_data(big_mask,:);

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

covariate_small_cg = cat(3,year_small_Thai);
covariate_small = cat(3,  madult_small_Thai, fadult_small_Thai, ...
                         kids_small_Thai, maleh_small_Thai, farm_small_Thai, ageh_small_Thai, age2h_small_Thai,educh_small_Thai);

covariate_big_cg = cat(3,year_big_Thai);
covariate_big = cat(3,  madult_big_Thai, fadult_big_Thai, ...
                         kids_big_Thai, maleh_big_Thai, farm_big_Thai, ageh_big_Thai, age2h_big_Thai, educh_big_Thai);

end

% Tuning Parameters
p = 1;

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
ngrid = 200;

% grid
grid = linspace(min(min(min(sample1),min(sample2)),[],'all'),max(max(max(sample1),max(sample2)),[],'all'),ngrid)';

% Estimation Starts
% Function for linear operator

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
    

for b = 1:200
disp(b)

% Bootstrapping 'Path-Wise'
index_for_path_wise_construnction1 = 0:N1:N1*T;
bootstrap_index1 = randi([1,N1],N1,1,btsp);
pw_btsp_index1 = bootstrap_index1 + repmat(index_for_path_wise_construnction1,1,1,btsp);

index_for_path_wise_construnction2 = 0:N2:N2*T;
bootstrap_index2 = randi([1,N2],N2,1,btsp);
pw_btsp_index2 = bootstrap_index2 + repmat(index_for_path_wise_construnction2,1,1,btsp);

% Check whether we have to conduct regression

if residual_flag == 1
    
id_index1 = sample1_id(bootstrap_index1,2);
id_index2 = sample2_id(bootstrap_index2,2);

id_index = [id_index1; id_index2];
id_index = unique(id_index);
index_mask = zeros(size(Thai_data.id));

for i = 1: size(id_index)
    temp = Thai_data.id == id_index(i);
    index_mask = index_mask + temp;
end
index_mask = logical(index_mask);

table_b = Thai_data(index_mask,:);

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

% Concatenating for Bootstrap Regression \hat \theta*

table_b1 = table();

table_b1.tc = y_b1;
table_b1.madult = covr_b1(:,1);                     
table_b1.fadult = covr_b1(:,2);                     
table_b1.kids = covr_b1(:,3);                     
table_b1.maleh = covr_b1(:,4);                     
table_b1.farm = covr_b1(:,5);                     
table_b1.ageh = covr_b1(:,6);                     
table_b1.age2h = covr_b1(:,7);                     
table_b1.year = covr_cg_b1(:,1);
table_b1.educh = covr_b1(:,8);


table_b2 = table();

table_b2.tc = y_b2;
table_b2.madult = covr_b2(:,1);                     
table_b2.fadult = covr_b2(:,2);                     
table_b2.kids = covr_b2(:,3);                     
table_b2.maleh = covr_b2(:,4);                     
table_b2.farm = covr_b2(:,5);                     
table_b2.ageh = covr_b2(:,6);                     
table_b2.age2h = covr_b2(:,7);                     
table_b2.year = covr_cg_b2(:,1);
table_b2.educh = covr_b2(:,8);

b_fitted = fitlm(table_b, variable_y + '~ year + madult + fadult + kids + maleh'...
                                + ' + farm + ageh + age2h + educh','Intercept',false);

b1sample1 = y_b1 - b_fitted.predict(table_b1);
b1sample1 = b1sample1(pw_btsp_index1);
b2sample2 = y_b2 - b_fitted.predict(table_b2);
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
