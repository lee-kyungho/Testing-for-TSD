clear
% TSD Test for evaluating the Carbon ETS in China.
% Kyungho Lee, Oliver Linton, Yoon-Jae Whang

% Read Data
NCDC_panel = readtable("NCDC_panel_quarterly.csv");
variables = ["visib" 'wdsp' "temp" "prcp" "provgb" "quarter"];
NCDC_panel = NCDC_panel(:,[variables, "year" 'CNT' 'stn']);

% Year Quarter
NCDC_panel.yq = NCDC_panel.year*10 + NCDC_panel.quarter;
NCDC_panel.time = (NCDC_panel.year - 2011)*4 + NCDC_panel.quarter;

Nan_index = logical(sum(isnan(table2array(NCDC_panel)),2));
NCDC_panel = NCDC_panel(~Nan_index,:);

% Dummy for Phase
NCDC_panel.phase0 = (NCDC_panel.year >= 2011) .* (NCDC_panel.year < 2014);
NCDC_panel.phase1 = (NCDC_panel.year >= 2014) .* (NCDC_panel.yq < 20163);
NCDC_panel.phase2 = (NCDC_panel.yq >= 20163);

% Flag for data allocation
NCDC_CNT0_ph0_flag = logical((NCDC_panel.CNT == 0) .* NCDC_panel.phase0 == 1);
NCDC_CNT0_ph1_flag = logical((NCDC_panel.CNT == 0) .* NCDC_panel.phase1 == 1);
NCDC_CNT0_ph2_flag = logical((NCDC_panel.CNT == 0) .* NCDC_panel.phase2 == 1);
NCDC_CNT1_ph0_flag = logical((NCDC_panel.CNT == 1) .* NCDC_panel.phase0 == 1);
NCDC_CNT1_ph1_flag = logical((NCDC_panel.CNT == 1) .* NCDC_panel.phase1 == 1);
NCDC_CNT1_ph2_flag = logical((NCDC_panel.CNT == 1) .* NCDC_panel.phase2 == 1);

% Allocating Data
NCDC_CNT0_ph0 = NCDC_panel(NCDC_CNT0_ph0_flag,:);
NCDC_CNT0_ph1 = NCDC_panel(NCDC_CNT0_ph1_flag,:);
NCDC_CNT0_ph2 = NCDC_panel(NCDC_CNT0_ph2_flag,:);
NCDC_CNT1_ph0 = NCDC_panel(NCDC_CNT1_ph0_flag,:);
NCDC_CNT1_ph1 = NCDC_panel(NCDC_CNT1_ph1_flag,:);
NCDC_CNT1_ph2 = NCDC_panel(NCDC_CNT1_ph2_flag,:);

% Time adjustmnet
NCDC_CNT0_ph2.time = NCDC_CNT0_ph2.time - max(NCDC_CNT0_ph1.time);
NCDC_CNT0_ph1.time = NCDC_CNT0_ph1.time - max(NCDC_CNT0_ph0.time);
NCDC_CNT1_ph2.time = NCDC_CNT1_ph2.time - max(NCDC_CNT1_ph1.time);
NCDC_CNT1_ph1.time = NCDC_CNT1_ph1.time - max(NCDC_CNT1_ph0.time);

% Concatenated Data
NCDC_ph0 = [NCDC_CNT0_ph0; NCDC_CNT1_ph0];
NCDC_ph1 = [NCDC_CNT0_ph1; NCDC_CNT1_ph1];
NCDC_ph2 = [NCDC_CNT0_ph2; NCDC_CNT1_ph2];
NCDC_phall = [NCDC_ph0; NCDC_ph1; NCDC_ph2];

NCDC_phall.provgb = categorical(NCDC_phall.provgb);
NCDC_phall.quarter = categorical(NCDC_phall.quarter);

% We use residuals for testing
residual_flag = 1;

% ------- TEST STARTS ---------- %

Results = ["design" 'residual_flag' "phase" 'order' 'contact'];

for phase = [0 1 2]

if  phase == 0
    sample_data  = NCDC_ph0;
    sample_data0 = NCDC_CNT0_ph0;
    sample_data1 = NCDC_CNT1_ph0;

elseif phase == 1
    sample_data  = NCDC_ph1;
    sample_data0 = NCDC_CNT0_ph1;
    sample_data1 = NCDC_CNT1_ph1;

elseif phase == 2
    sample_data  = NCDC_ph2;
    sample_data0 = NCDC_CNT0_ph2;
    sample_data1 = NCDC_CNT1_ph2;
end

sample_data.provgb = categorical(sample_data.provgb);
sample_data.quarter = categorical(sample_data.quarter);

[visib_test_data_0, unique_id_0] = test_data(sample_data0, "visib");
[visib_test_data_1, unique_id_1] = test_data(sample_data1, "visib");

ngrid = 100;

% Tuning Parameters
p = 1;

% Bootstrapping
btsp = 1;

if residual_flag == 1

% Linear Regression
fitted = fitlm(NCDC_phall, "visib" + '~ wdsp + temp + prcp + provgb + quarter'); % OLS regression

res_visib_data_0 = sample_data0.visib - fitted.predict(sample_data(sample_data.CNT == 0, :));
res_visib_data_1 = sample_data1.visib - fitted.predict(sample_data(sample_data.CNT == 1, :));

sample_data0.res = res_visib_data_0;
sample_data1.res = res_visib_data_1;

[res_visib_test_data_0, ~] = test_data(sample_data0, "res");
[res_visib_test_data_1, ~] = test_data(sample_data1, "res");

[wdsp_test_data_0, ~] = test_data(sample_data0, "wdsp");
[wdsp_test_data_1, ~] = test_data(sample_data1, "wdsp");

[temp_test_data_0, ~] = test_data(sample_data0, "temp");
[temp_test_data_1, ~] = test_data(sample_data1, "temp");

[prcp_test_data_0, ~] = test_data(sample_data0, "prcp");
[prcp_test_data_1, ~] = test_data(sample_data1, "prcp");

[provgb_test_data_0, ~] = test_data(sample_data0, "provgb");
[provgb_test_data_1, ~] = test_data(sample_data1, "provgb");

[quarter_test_data_0, ~] = test_data(sample_data0, "quarter");
[quarter_test_data_1, ~] = test_data(sample_data1, "quarter");

covariate_0 = cat(3,  wdsp_test_data_0, temp_test_data_0, ...
                         prcp_test_data_0, provgb_test_data_0, quarter_test_data_0);

covariate_1 = cat(3,  wdsp_test_data_1, temp_test_data_1, ...
                         prcp_test_data_1, provgb_test_data_1, quarter_test_data_1);
end

for design = [1 2]

    
disp(design)    

if design == 1
    disp("Not Regulated >= Regulated")    
    
if residual_flag == 1
    sample1 = squeeze(res_visib_test_data_0(:,:));
    sample2 = squeeze(res_visib_test_data_1(:,:));
    sample1_id = unique_id_0;
    sample2_id = unique_id_1;
    
    sample1_covariate = covariate_0;
    sample2_covariate = covariate_1;

else
    sample1 = squeeze(visib_test_data_0(:,:));
    sample2 = squeeze(visib_test_data_1(:,:));

end

elseif design == 2
    
    disp("Regulated >= Not Regulated")    
if residual_flag == 1
    sample1 = squeeze(res_visib_test_data_1(:,:));
    sample2 = squeeze(res_visib_test_data_0(:,:));
    sample1_id = unique_id_1;
    sample2_id = unique_id_0;
    
    sample1_covariate = covariate_1;
    sample2_covariate = covariate_0;

else
    sample1 = squeeze(visib_test_data_1(:,:));
    sample2 = squeeze(visib_test_data_0(:,:));
end

end

ngrid = 100;
% grid
grid = linspace(min(min(min(sample1),min(sample2)),[],'all'),max(max(max(sample1),max(sample2)),[],'all'),ngrid)';

% Estimation Starts
T = size(sample1,2)-1;
N1 = size(sample1,1);
N2 = size(sample2,1);

r_N = sqrt(N1*N2/(N1+N2));

T_11_max_contact_2_list = [];
T_12_max_contact_2_list = [];
T_21_max_contact_2_list = [];
T_22_max_contact_2_list = [];

% (Time, SD) order = (1,1)
op1_11 = operation(1,1,sample1,grid); 
op2_11 = operation(1,1,sample2,grid); 
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

Lamb_11_max = Lambda(D_11,p,'max'); % Output dim: ngrid * 1
Lamb_12_max = Lambda(D_12,p,'max'); % Output dim: ngrid * 1
Lamb_21_max = Lambda(D_21,p,'max'); % Output dim: ngrid * 1

% For n=2, m=2, we need to make collection of D_22 and D_12_T
Lamb_22_max = Lambda(D_22_collection,p,'max'); % Output dim: ngrid * 1

%---- Integration for T_{N} ----------%

T_11_max = r_N^p * trapz(Lamb_11_max);
T_12_max = r_N^p * trapz(Lamb_12_max);
T_21_max = r_N^p * trapz(Lamb_21_max);
T_22_max = r_N^p * trapz(Lamb_22_max);


parfor b = 1:200
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
    
id_index1 = sample1_id(bootstrap_index1);
id_index2 = sample2_id(bootstrap_index2);

id_index = [id_index1; id_index2];
%id_index = unique(id_index);
index_mask = zeros(size(NCDC_phall.stn));

for i = 1: size(id_index,1)
    temp = NCDC_phall.stn == id_index(i);
    index_mask = index_mask + temp;
end

index_logic = logical(index_mask);
index_wgt   = index_mask(index_logic ~= 0);

table_b = NCDC_phall(index_logic,variables);
table_b.provgb = categorical(table_b.provgb);
table_b.quarter = categorical(table_b.quarter);

% sample 1
y_bsample1 = sample1(pw_btsp_index1);
covr_bsample1 = sample1_covariate(bootstrap_index1,:,:);

% Unstacking for regression
y_b1 = reshape(y_bsample1,[size(y_bsample1,1)*size(y_bsample1,2),1]);
covr_b1 = reshape(covr_bsample1,[size(covr_bsample1,1)*size(covr_bsample1,2),size(covr_bsample1,3)]);

% sample 2
y_bsample2 = sample2(pw_btsp_index2);
covr_bsample2 = sample2_covariate(bootstrap_index2,:,:);

% Unstacking for regression
y_b2 = reshape(y_bsample2,[size(y_bsample2,1)*size(y_bsample2,2),1]);
covr_b2 = reshape(covr_bsample2,[size(covr_bsample2,1)*size(covr_bsample2,2),size(covr_bsample2,3)]);

% Concatenating for Bootstrap Regression \hat \theta*
table_b1 = table();

table_b1.visib = y_b1;
table_b1.wdsp = covr_b1(:,1);                     
table_b1.temp = covr_b1(:,2);                     
table_b1.prcp = covr_b1(:,3);                     
table_b1.provgb = categorical(covr_b1(:,4));                     
table_b1.quarter = categorical(covr_b1(:,5));  

table_b2 = table();

table_b2.visib = y_b2;
table_b2.wdsp = covr_b2(:,1);                     
table_b2.temp = covr_b2(:,2);                     
table_b2.prcp = covr_b2(:,3);                     
table_b2.provgb = categorical(covr_b2(:,4));                     
table_b2.quarter = categorical(covr_b2(:,5));                     

% Note on Warning
% Depending on bootstrap draws, some categorical values are multicollinear
% One dummy will be dropped while fitlm runs
b_fitted = fitlm(table_b, "visib" + '~ wdsp + temp + prcp + provgb + quarter','Weight', index_wgt);

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

% For n = 2, m = 2, we need to make collection of D_22 and D_12_T
% Reshaping for concatenation
b_D_12_T_recentered = reshape(b_D_12_T_recentered,[size(b_D_12_T_recentered,1),1,size(b_D_12_T_recentered,2)]);
b_D_22_collection_recentered = cat(2,b_D_22_recentered,b_D_12_T_recentered); % ngrid * J * b

%---- Bootsrap Sample Statistics -----%
%---- Contact Set Approach -----------%

% Setting Tuning parameter 
N = (N1 + N2) / 2;

% c = 0.7
c_N_2 = 0.7*log(log(N));

[T_11_max_contact_2, T_12_max_contact_2, ...
 T_21_max_contact_2, T_22_max_contact_2] = ...
            TSD_contact_China(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p, N, r_N, c_N_2);


T_11_max_contact_2_list = [T_11_max_contact_2_list T_11_max_contact_2];
T_12_max_contact_2_list = [T_12_max_contact_2_list T_12_max_contact_2];
T_21_max_contact_2_list = [T_21_max_contact_2_list T_21_max_contact_2];
T_22_max_contact_2_list = [T_22_max_contact_2_list T_22_max_contact_2];

end

% % Contact Set Approach
p_value_11_max_contact_2 = mean(T_11_max <= T_11_max_contact_2_list);
p_value_12_max_contact_2 = mean(T_12_max <= T_12_max_contact_2_list);
p_value_21_max_contact_2 = mean(T_21_max <= T_21_max_contact_2_list);
p_value_22_max_contact_2 = mean(T_22_max <= T_22_max_contact_2_list);

Results = [Results; [design residual_flag phase "(1,1)" p_value_11_max_contact_2;
                     design residual_flag phase "(1,2)" p_value_12_max_contact_2;
                     design residual_flag phase "(2,1)" p_value_21_max_contact_2;
                     design residual_flag phase "(2,2)" p_value_22_max_contact_2;
                    ]]

end
writematrix(Results, "China_CNT_1028.csv")
end

