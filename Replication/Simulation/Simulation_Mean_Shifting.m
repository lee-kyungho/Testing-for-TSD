clear

% Simulation
% Author: Kyungho Lee, Oliver Linton, and Yoon-Jae Whang
% Interior Distributions Size / Power
% Mean-Shifitng DGP of Uniform Distributions
% We change the degree of mean shifted 

T = 4; % Terminal Periods
ngrid = 100; % Number of Grid Points
btsp = 200; % Number of bootstrapiing
n_simul = 1000; % the number of simulation
p1=1;  %  L1 staistics
p2=2;  %  L2 staistics

% Significance Level
alpha = 0.05;

% Tuning param eta
eta = 10^(-6);

Results_max_p1 = ["N" "design" "order" "LFC" ...
                   "contact_1" "contact_2" "contact_3" "contact_4" "contact_5"  ... 
                   "NDM_1" "NDM_2" "NDM_3" "NDM_4" "NDM_5" ...
                    ];

Results_sum_p1 = ["N" "design" "order" "LFC" ...
                   "contact_1" "contact_2" "contact_3" "contact_4" "contact_5"  ... 
                   "NDM_1" "NDM_2" "NDM_3" "NDM_4" "NDM_5" ...
                   ];

Results_max_p2 = [ "N" "design" "(2,2)" "LFC" ...
                   "contact_1" "contact_2" "contact_3" "contact_4" "contact_5"... 
                   "NDM_1_1" "NDM_1_2" "NDM_1_3" "NDM_1_4" "NDM_1_5" ...
                   "NDM_2_1" "NDM_2_2" "NDM_2_3" "NDM_2_4" "NDM_2_5" ];

              
Results_sum_p2 =  [ "N" "design" "(2,2)" "LFC" ...
                   "contact_1" "contact_2" "contact_3" "contact_4" "contact_5" ... 
                   "NDM_1_1" "NDM_1_2" "NDM_1_3" "NDM_1_4" "NDM_1_5" ...
                   "NDM_2_1" "NDM_2_2" "NDM_2_3" "NDM_2_4" "NDM_2_5" 
                   ];

% Run Simulation
for N = [100, 500] 
for design = 0:30
    
%---------%--------- DGP --------%-----------%--------%
X = zeros(N,T+1,n_simul);
X0 = unifrnd(0,1,N,1,n_simul);
X(:,1,:) = X0;

for i = 1:T
    X(:,i+1,:) = 0.5*X(:,i,:) +  unifrnd(0,1,N,1,n_simul)  - 1/2;
end

Y = zeros(N,T+1,n_simul);
Y0 = unifrnd(0,1,N,1,n_simul);
Y(:,1,:) = Y0;

for i = 1:T
    Y(:,i+1,:) = 0.5*Y(:,i,:) + unifrnd(0,1,N,1,n_simul) - 1/2;
end

% Shfting distributions
% design 0 corresponds to the least favorable case

if design == 1
Y(:,1,:) = Y(:,1,:) + 0.02;
elseif design == 2
Y(:,1,:) = Y(:,1,:) + 0.04;
elseif design == 3
Y(:,1,:) = Y(:,1,:) + 0.06;
elseif design == 4
Y(:,1,:) = Y(:,1,:) + 0.08;
elseif design == 5
Y(:,1,:) = Y(:,1,:) + 0.10;
elseif design == 6
Y(:,1,:) = Y(:,1,:) + 0.12;
elseif design == 7
Y(:,1,:) = Y(:,1,:) + 0.14;
elseif design == 8
Y(:,1,:) = Y(:,1,:) - 0.02;
elseif design == 9
Y(:,1,:) = Y(:,1,:) - 0.06;
elseif design == 10
Y(:,1,:) = Y(:,1,:) - 0.10;
elseif design == 11
Y(:,3,:) = Y(:,3,:) + 0.02;
elseif design == 12
Y(:,3,:) = Y(:,3,:) + 0.04;
elseif design == 13
Y(:,3,:) = Y(:,3,:) + 0.06;
elseif design == 14
Y(:,3,:) = Y(:,3,:) + 0.08;
elseif design == 15
Y(:,3,:) = Y(:,3,:) + 0.10;
elseif design == 16
Y(:,3,:) = Y(:,3,:) + 0.12;
elseif design == 17
Y(:,3,:) = Y(:,3,:) + 0.14;
elseif design == 18
Y(:,3,:) = Y(:,3,:) - 0.02;
elseif design == 19
Y(:,3,:) = Y(:,3,:) - 0.06;
elseif design == 20
Y(:,3,:) = Y(:,3,:) - 0.10;
elseif design == 21
Y(:,T+1,:) = Y(:,T+1,:) + 0.02;
elseif design == 22
Y(:,T+1,:) = Y(:,T+1,:) + 0.04;
elseif design == 23
Y(:,T+1,:) = Y(:,T+1,:) + 0.06;
elseif design == 24
Y(:,T+1,:) = Y(:,T+1,:) + 0.08;
elseif design == 25
Y(:,T+1,:) = Y(:,T+1,:) + 0.10;
elseif design == 26
Y(:,T+1,:) = Y(:,T+1,:) + 0.12;
elseif design == 27
Y(:,T+1,:) = Y(:,T+1,:) + 0.14;
elseif design == 28
Y(:,T+1,:) = Y(:,T+1,:) - 0.02;
elseif design == 29
Y(:,T+1,:) = Y(:,T+1,:) - 0.06;
elseif design == 30
Y(:,T+1,:) = Y(:,T+1,:) - 0.10;
end
% for two differenet sample size: 
% r_N = sqrt(size(X,1)*size(Y,1)/(size(X,1) + size(Y,1)));
r_N = sqrt(size(X,1));

%---------%---------%---------%---------%---------%---------%---------
%---------%---------% Simulation Starts %---------%---------%---------
%---------%---------%---------%---------%---------%---------%---------

% alpha = 0.05
%(1,1)

% p = 1
% p = 1 LFC

rejection_11_btsp_max_LFC_p1     = []; 
rejection_11_btsp_sum_LFC_p1     = [];

% p = 1 Contact

rejection_11_max_contact_1_p1 = [];
rejection_11_sum_contact_1_p1 = [];

rejection_11_max_contact_2_p1 = [];
rejection_11_sum_contact_2_p1 = [];

rejection_11_max_contact_3_p1 = [];
rejection_11_sum_contact_3_p1 = [];

rejection_11_max_contact_4_p1 = [];
rejection_11_sum_contact_4_p1 = [];

rejection_11_max_contact_5_p1 = [];
rejection_11_sum_contact_5_p1 = [];

% p = 1 NDM

rejection_11_btsp_max_NDM_1_p1 = [];
rejection_11_btsp_sum_NDM_1_p1 = [];

rejection_11_btsp_max_NDM_2_p1 = [];
rejection_11_btsp_sum_NDM_2_p1 = [];

rejection_11_btsp_max_NDM_3_p1 = [];
rejection_11_btsp_sum_NDM_3_p1 = [];

rejection_11_btsp_max_NDM_4_p1 = [];
rejection_11_btsp_sum_NDM_4_p1 = [];

rejection_11_btsp_max_NDM_5_p1 = [];
rejection_11_btsp_sum_NDM_5_p1 = [];

% p = 2 contact

rejection_11_btsp_max_LFC_p2     = []; 
rejection_11_btsp_sum_LFC_p2     = [];

rejection_11_max_contact_1_p2 = [];
rejection_11_sum_contact_1_p2 = [];

rejection_11_max_contact_2_p2 = [];
rejection_11_sum_contact_2_p2 = [];

rejection_11_max_contact_3_p2 = [];
rejection_11_sum_contact_3_p2 = [];

rejection_11_max_contact_4_p2 = [];
rejection_11_sum_contact_4_p2 = [];

rejection_11_max_contact_5_p2 = [];
rejection_11_sum_contact_5_p2 = [];

% p = 2 First way of Approximation

rejection_11_btsp_max_NDM_1_p2_first = [];
rejection_11_btsp_sum_NDM_1_p2_first = [];

rejection_11_btsp_max_NDM_2_p2_first = [];
rejection_11_btsp_sum_NDM_2_p2_first = [];

rejection_11_btsp_max_NDM_3_p2_first = [];
rejection_11_btsp_sum_NDM_3_p2_first = [];

rejection_11_btsp_max_NDM_4_p2_first = [];
rejection_11_btsp_sum_NDM_4_p2_first = [];

rejection_11_btsp_max_NDM_5_p2_first = [];
rejection_11_btsp_sum_NDM_5_p2_first = [];

% p = 2 Second way of Approximation

rejection_11_btsp_max_NDM_1_p2_second = [];
rejection_11_btsp_sum_NDM_1_p2_second = [];

rejection_11_btsp_max_NDM_2_p2_second = [];
rejection_11_btsp_sum_NDM_2_p2_second = [];

rejection_11_btsp_max_NDM_3_p2_second = [];
rejection_11_btsp_sum_NDM_3_p2_second = [];

rejection_11_btsp_max_NDM_4_p2_second = [];
rejection_11_btsp_sum_NDM_4_p2_second = [];

rejection_11_btsp_max_NDM_5_p2_second = [];
rejection_11_btsp_sum_NDM_5_p2_second = [];

% (1,2)

% p = 1 LFC

rejection_12_btsp_max_LFC_p1     = []; 
rejection_12_btsp_sum_LFC_p1     = [];

% p = 1 Contact

rejection_12_max_contact_1_p1 = [];
rejection_12_sum_contact_1_p1 = [];

rejection_12_max_contact_2_p1 = [];
rejection_12_sum_contact_2_p1 = [];

rejection_12_max_contact_3_p1 = [];
rejection_12_sum_contact_3_p1 = [];

rejection_12_max_contact_4_p1 = [];
rejection_12_sum_contact_4_p1 = [];

rejection_12_max_contact_5_p1 = [];
rejection_12_sum_contact_5_p1 = [];

% p = 1 NDM 

rejection_12_btsp_max_NDM_1_p1 = [];
rejection_12_btsp_sum_NDM_1_p1 = [];

rejection_12_btsp_max_NDM_2_p1 = [];
rejection_12_btsp_sum_NDM_2_p1 = [];

rejection_12_btsp_max_NDM_3_p1 = [];
rejection_12_btsp_sum_NDM_3_p1 = [];

rejection_12_btsp_max_NDM_4_p1 = [];
rejection_12_btsp_sum_NDM_4_p1 = [];

rejection_12_btsp_max_NDM_5_p1 = [];
rejection_12_btsp_sum_NDM_5_p1 = [];

% p = 2 LFC

rejection_12_btsp_max_LFC_p2     = []; 
rejection_12_btsp_sum_LFC_p2     = [];

% p = 2 Contact

rejection_12_max_contact_1_p2 = [];
rejection_12_sum_contact_1_p2 = [];

rejection_12_max_contact_2_p2 = [];
rejection_12_sum_contact_2_p2 = [];

rejection_12_max_contact_3_p2 = [];
rejection_12_sum_contact_3_p2 = [];

rejection_12_max_contact_4_p2 = [];
rejection_12_sum_contact_4_p2 = [];

rejection_12_max_contact_5_p2 = [];
rejection_12_sum_contact_5_p2 = [];

% p = 2 NDM First

rejection_12_btsp_max_NDM_1_p2_first = [];
rejection_12_btsp_sum_NDM_1_p2_first = [];

rejection_12_btsp_max_NDM_2_p2_first = [];
rejection_12_btsp_sum_NDM_2_p2_first = [];

rejection_12_btsp_max_NDM_3_p2_first = [];
rejection_12_btsp_sum_NDM_3_p2_first = [];

rejection_12_btsp_max_NDM_4_p2_first = [];
rejection_12_btsp_sum_NDM_4_p2_first = [];

rejection_12_btsp_max_NDM_5_p2_first = [];
rejection_12_btsp_sum_NDM_5_p2_first = [];

% p = 2 NDM Second

rejection_12_btsp_max_NDM_1_p2_second = [];
rejection_12_btsp_sum_NDM_1_p2_second = [];

rejection_12_btsp_max_NDM_2_p2_second = [];
rejection_12_btsp_sum_NDM_2_p2_second = [];

rejection_12_btsp_max_NDM_3_p2_second = [];
rejection_12_btsp_sum_NDM_3_p2_second = [];

rejection_12_btsp_max_NDM_4_p2_second = [];
rejection_12_btsp_sum_NDM_4_p2_second = [];

rejection_12_btsp_max_NDM_5_p2_second = [];
rejection_12_btsp_sum_NDM_5_p2_second = [];

% (2,1)

% p = 1 LFC

rejection_21_btsp_max_LFC_p1     = []; 
rejection_21_btsp_sum_LFC_p1     = [];

% p = 1 Contact

rejection_21_max_contact_1_p1 = [];
rejection_21_sum_contact_1_p1 = [];

rejection_21_max_contact_2_p1 = [];
rejection_21_sum_contact_2_p1 = [];

rejection_21_max_contact_3_p1 = [];
rejection_21_sum_contact_3_p1 = [];

rejection_21_max_contact_4_p1 = [];
rejection_21_sum_contact_4_p1 = [];

rejection_21_max_contact_5_p1 = [];
rejection_21_sum_contact_5_p1 = [];

% p = 1 NDM

rejection_21_btsp_max_NDM_1_p1 = [];
rejection_21_btsp_sum_NDM_1_p1 = [];

rejection_21_btsp_max_NDM_2_p1 = [];
rejection_21_btsp_sum_NDM_2_p1 = [];

rejection_21_btsp_max_NDM_3_p1 = [];
rejection_21_btsp_sum_NDM_3_p1 = [];

rejection_21_btsp_max_NDM_4_p1 = [];
rejection_21_btsp_sum_NDM_4_p1 = [];

rejection_21_btsp_max_NDM_5_p1 = [];
rejection_21_btsp_sum_NDM_5_p1 = [];
 
% p = 2 LFC

rejection_21_btsp_max_LFC_p2     = []; 
rejection_21_btsp_sum_LFC_p2     = [];

% p = 2 Contact

rejection_21_max_contact_1_p2 = [];
rejection_21_sum_contact_1_p2 = [];

rejection_21_max_contact_2_p2 = [];
rejection_21_sum_contact_2_p2 = [];

rejection_21_max_contact_3_p2 = [];
rejection_21_sum_contact_3_p2 = [];

rejection_21_max_contact_4_p2 = [];
rejection_21_sum_contact_4_p2 = [];

rejection_21_max_contact_5_p2 = [];
rejection_21_sum_contact_5_p2 = [];

% p = 2 NDM First

rejection_21_btsp_max_NDM_1_p2_first = [];
rejection_21_btsp_sum_NDM_1_p2_first = [];

rejection_21_btsp_max_NDM_2_p2_first = [];
rejection_21_btsp_sum_NDM_2_p2_first = [];

rejection_21_btsp_max_NDM_3_p2_first = [];
rejection_21_btsp_sum_NDM_3_p2_first = [];

rejection_21_btsp_max_NDM_4_p2_first = [];
rejection_21_btsp_sum_NDM_4_p2_first = [];

rejection_21_btsp_max_NDM_5_p2_first = [];
rejection_21_btsp_sum_NDM_5_p2_first = [];

% p = 2 NDM Second

rejection_21_btsp_max_NDM_1_p2_second = [];
rejection_21_btsp_sum_NDM_1_p2_second = [];

rejection_21_btsp_max_NDM_2_p2_second = [];
rejection_21_btsp_sum_NDM_2_p2_second = [];

rejection_21_btsp_max_NDM_3_p2_second = [];
rejection_21_btsp_sum_NDM_3_p2_second = [];

rejection_21_btsp_max_NDM_4_p2_second = [];
rejection_21_btsp_sum_NDM_4_p2_second = [];

rejection_21_btsp_max_NDM_5_p2_second = [];
rejection_21_btsp_sum_NDM_5_p2_second = [];

% (2,2)

% p = 1 LFC

rejection_22_btsp_max_LFC_p1     = []; 
rejection_22_btsp_sum_LFC_p1     = [];

% p = 1 Contact

rejection_22_max_contact_1_p1 = [];
rejection_22_sum_contact_1_p1 = [];

rejection_22_max_contact_2_p1 = [];
rejection_22_sum_contact_2_p1 = [];

rejection_22_max_contact_3_p1 = [];
rejection_22_sum_contact_3_p1 = [];

rejection_22_max_contact_4_p1 = [];
rejection_22_sum_contact_4_p1 = [];

rejection_22_max_contact_5_p1 = [];
rejection_22_sum_contact_5_p1 = [];

% p = 1 NDM 

rejection_22_btsp_max_NDM_1_p1 = [];
rejection_22_btsp_sum_NDM_1_p1 = [];

rejection_22_btsp_max_NDM_2_p1 = [];
rejection_22_btsp_sum_NDM_2_p1 = [];

rejection_22_btsp_max_NDM_3_p1 = [];
rejection_22_btsp_sum_NDM_3_p1 = [];

rejection_22_btsp_max_NDM_4_p1 = [];
rejection_22_btsp_sum_NDM_4_p1 = [];

rejection_22_btsp_max_NDM_5_p1 = [];
rejection_22_btsp_sum_NDM_5_p1 = [];

% p = 2 LFC

rejection_22_btsp_max_LFC_p2     = []; 
rejection_22_btsp_sum_LFC_p2     = [];

% p = 2 Contact

rejection_22_max_contact_1_p2 = [];
rejection_22_sum_contact_1_p2 = [];

rejection_22_max_contact_2_p2 = [];
rejection_22_sum_contact_2_p2 = [];

rejection_22_max_contact_3_p2 = [];
rejection_22_sum_contact_3_p2 = [];

rejection_22_max_contact_4_p2 = [];
rejection_22_sum_contact_4_p2 = [];

rejection_22_max_contact_5_p2 = [];
rejection_22_sum_contact_5_p2 = [];

% p = 2 NDM First

rejection_22_btsp_max_NDM_1_p2_first = [];
rejection_22_btsp_sum_NDM_1_p2_first = [];

rejection_22_btsp_max_NDM_2_p2_first = [];
rejection_22_btsp_sum_NDM_2_p2_first = [];

rejection_22_btsp_max_NDM_3_p2_first = [];
rejection_22_btsp_sum_NDM_3_p2_first = [];

rejection_22_btsp_max_NDM_4_p2_first = [];
rejection_22_btsp_sum_NDM_4_p2_first = [];

rejection_22_btsp_max_NDM_5_p2_first = [];
rejection_22_btsp_sum_NDM_5_p2_first = [];

% p = 2 NDM Second

rejection_22_btsp_max_NDM_1_p2_second = [];
rejection_22_btsp_sum_NDM_1_p2_second = [];

rejection_22_btsp_max_NDM_2_p2_second = [];
rejection_22_btsp_sum_NDM_2_p2_second = [];

rejection_22_btsp_max_NDM_3_p2_second = [];
rejection_22_btsp_sum_NDM_3_p2_second = [];

rejection_22_btsp_max_NDM_4_p2_second = [];
rejection_22_btsp_sum_NDM_4_p2_second = [];

rejection_22_btsp_max_NDM_5_p2_second = [];
rejection_22_btsp_sum_NDM_5_p2_second = [];

tic

disp("Simulation Starts")

parfor R = 1:n_simul

disp("Simulation:")
disp(R)
    
    
sample1 = squeeze(X(:,:,R));
sample2 = squeeze(Y(:,:,R));


% grid
grid = linspace(min(min(sample1,sample2),[],'all'),max(max(sample1,sample2),[],'all'),ngrid)';
% Function for linear operator

% Bootstrapping 'Path-Wise'
b1sample1 = path_wise_bootstrap(sample1,btsp);
b2sample2 = path_wise_bootstrap(sample2,btsp);

% Calculate Test Statistics
% Note on dimension on D

% (Time, SD) order = (1,1)
op1_11 = operation(1,1,sample1,grid); % Output dim ngrid * (T+1) * btsp
op2_11 = operation(1,1,sample2,grid); % Output dim ngrid * (T+1) * btsp
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

% p = 1 (p1)

Lamb_11_max_p1 = Lambda(D_11,p1,'max'); % Output dim: ngrid * 1
Lamb_11_sum_p1 = Lambda(D_11,p1,'sum'); % Output dim: ngrid * 1

Lamb_12_max_p1 = Lambda(D_12,p1,'max'); % Output dim: ngrid * 1
Lamb_12_sum_p1 = Lambda(D_12,p1,'sum'); % Output dim: ngrid * 1

Lamb_21_max_p1 = Lambda(D_21,p1,'max'); % Output dim: ngrid * 1
Lamb_21_sum_p1 = Lambda(D_21,p1,'sum'); % Output dim: ngrid * 1

% For n=2, m=2, we need to make collection of D_22 and D_12_T

D_22_collection = cat(2,D_22,D_12_T); % ngrid * J

Lamb_22_max_p1 = Lambda(D_22_collection,p1,'max'); % Output dim: ngrid * 1
Lamb_22_sum_p1 = Lambda(D_22_collection,p1,'sum'); % Output dim: ngrid * 1

% p = 2 (p2)

Lamb_11_max_p2 = Lambda(D_11,p2,'max'); % Output dim: ngrid * 1
Lamb_11_sum_p2 = Lambda(D_11,p2,'sum'); % Output dim: ngrid * 1

Lamb_12_max_p2 = Lambda(D_12,p2,'max'); % Output dim: ngrid * 1
Lamb_12_sum_p2 = Lambda(D_12,p2,'sum'); % Output dim: ngrid * 1

Lamb_21_max_p2 = Lambda(D_21,p2,'max'); % Output dim: ngrid * 1
Lamb_21_sum_p2 = Lambda(D_21,p2,'sum'); % Output dim: ngrid * 1

Lamb_22_max_p2 = Lambda(D_22_collection,p2,'max'); % Output dim: ngrid * 1
Lamb_22_sum_p2 = Lambda(D_22_collection,p2,'sum'); % Output dim: ngrid * 1

%---- sumegration for T_{N} ----------%

T_11_max_p1 = r_N^p1 * trapz(Lamb_11_max_p1);
T_11_sum_p1 = r_N^p1 * trapz(Lamb_11_sum_p1);

T_12_max_p1 = r_N^p1 * trapz(Lamb_12_max_p1);
T_12_sum_p1 = r_N^p1 * trapz(Lamb_12_sum_p1);

T_21_max_p1 = r_N^p1 * trapz(Lamb_21_max_p1);
T_21_sum_p1 = r_N^p1 * trapz(Lamb_21_sum_p1);

T_22_max_p1 = r_N^p1 * trapz(Lamb_22_max_p1);
T_22_sum_p1 = r_N^p1 * trapz(Lamb_22_sum_p1);

% p = 2

T_11_max_p2 = r_N^p2 * trapz(Lamb_11_max_p2);
T_11_sum_p2 = r_N^p2 * trapz(Lamb_11_sum_p2);

T_12_max_p2 = r_N^p2 * trapz(Lamb_12_max_p2);
T_12_sum_p2 = r_N^p2 * trapz(Lamb_12_sum_p2);

T_21_max_p2 = r_N^p2 * trapz(Lamb_21_max_p2);
T_21_sum_p2 = r_N^p2 * trapz(Lamb_21_sum_p2);

T_22_max_p2 = r_N^p2 * trapz(Lamb_22_max_p2);
T_22_sum_p2 = r_N^p2 * trapz(Lamb_22_sum_p2);


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

% Bootstrap Sample ------------------%

% p = 1

btsp_Lamb_11_max_LFC_p1 = Lambda(b_D_11_recentered,p1,'max'); % Output dim: ngrid * 1 * btsp
btsp_Lamb_11_sum_LFC_p1 = Lambda(b_D_11_recentered,p1,'sum'); % Output dim: ngrid * 1 * btsp

btsp_Lamb_12_max_LFC_p1 = Lambda(b_D_12_recentered,p1,'max'); % Output dim: ngrid * 1 * btsp
btsp_Lamb_12_sum_LFC_p1 = Lambda(b_D_12_recentered,p1,'sum'); % Output dim: ngrid * 1 * btsp

btsp_Lamb_21_max_LFC_p1 = Lambda(b_D_21_recentered,p1,'max'); % Output dim: ngrid * 1 * btsp
btsp_Lamb_21_sum_LFC_p1 = Lambda(b_D_21_recentered,p1,'sum'); % Output dim: ngrid * 1 * btsp

% For n = 2, m = 2, we need to make collection of D_22 and D_12_T
% Reshaping for concatenation
b_D_12_T_recentered = reshape(b_D_12_T_recentered,[size(b_D_12_T_recentered,1),1,size(b_D_12_T_recentered,2)]);

b_D_22_collection_recentered = cat(2,b_D_22_recentered,b_D_12_T_recentered); % ngrid * J * btsp

btsp_Lamb_22_max_LFC_p1 = Lambda(b_D_22_collection_recentered,p1,'max'); % Output dim: ngrid * 1 * btsp
btsp_Lamb_22_sum_LFC_p1 = Lambda(b_D_22_collection_recentered,p1,'sum'); % Output dim: ngrid * 1 * btsp

% p = 2

btsp_Lamb_11_max_LFC_p2 = Lambda(b_D_11_recentered,p2,'max'); % Output dim: ngrid * 1 * btsp
btsp_Lamb_11_sum_LFC_p2 = Lambda(b_D_11_recentered,p2,'sum'); % Output dim: ngrid * 1 * btsp

btsp_Lamb_12_max_LFC_p2 = Lambda(b_D_12_recentered,p2,'max'); % Output dim: ngrid * 1 * btsp
btsp_Lamb_12_sum_LFC_p2 = Lambda(b_D_12_recentered,p2,'sum'); % Output dim: ngrid * 1 * btsp

btsp_Lamb_21_max_LFC_p2 = Lambda(b_D_21_recentered,p2,'max'); % Output dim: ngrid * 1 * btsp
btsp_Lamb_21_sum_LFC_p2 = Lambda(b_D_21_recentered,p2,'sum'); % Output dim: ngrid * 1 * btsp

btsp_Lamb_22_max_LFC_p2 = Lambda(b_D_22_collection_recentered,p2,'max'); % Output dim: ngrid * 1 * btsp
btsp_Lamb_22_sum_LFC_p2 = Lambda(b_D_22_collection_recentered,p2,'sum'); % Output dim: ngrid * 1 * btsp

%---- Bootstrap version Integral of T_{N} %

% Bootstrap Sample
% p = 1

btsp_T_11_max_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_11_max_LFC_p1,1);
btsp_T_11_sum_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_11_sum_LFC_p1,1);

btsp_T_12_max_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_12_max_LFC_p1,1);
btsp_T_12_sum_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_12_sum_LFC_p1,1);

btsp_T_21_max_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_21_max_LFC_p1,1);
btsp_T_21_sum_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_21_sum_LFC_p1,1);

btsp_T_22_max_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_22_max_LFC_p1,1);
btsp_T_22_sum_LFC_p1 = r_N^p1 * trapz(btsp_Lamb_22_sum_LFC_p1,1);

% p = 2

btsp_T_11_max_LFC_p2 = r_N^p2 * trapz(btsp_Lamb_11_max_LFC_p2,1);
btsp_T_11_sum_LFC_p2 = r_N^p2 * trapz(btsp_Lamb_11_sum_LFC_p2,1);

btsp_T_12_max_LFC_p2 = r_N^p2 * trapz(btsp_Lamb_12_max_LFC_p2,1);
btsp_T_12_sum_LFC_p2 = r_N^p2 * trapz(btsp_Lamb_12_sum_LFC_p2,1);

btsp_T_21_max_LFC_p2 = r_N^p2 * trapz(btsp_Lamb_21_max_LFC_p2,1);
btsp_T_21_sum_LFC_p2 = r_N^p2 * trapz(btsp_Lamb_21_sum_LFC_p2,1);

btsp_T_22_max_LFC_p2 = r_N^p2 * trapz(btsp_Lamb_22_max_LFC_p2,1);
btsp_T_22_sum_LFC_p2 = r_N^p2 * trapz(btsp_Lamb_22_sum_LFC_p2,1);

%---- Bootsrap Sample Statistics -----%
%---- Contact Set Approach -----------%

% Setting Tuning parameter 
% Tuning parameters range 0.3 ~ 0.6
% 1. c = 0.3
c_N_1 = 0.3*log(log(N));

[T_11_max_contact_1_p1, T_11_sum_contact_1_p1, T_12_max_contact_1_p1, T_12_sum_contact_1_p1, ...
 T_21_max_contact_1_p1, T_21_sum_contact_1_p1, T_22_max_contact_1_p1, T_22_sum_contact_1_p1] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p1, N, r_N, c_N_1);

                       
[T_11_max_contact_1_p2, T_11_sum_contact_1_p2, T_12_max_contact_1_p2, T_12_sum_contact_1_p2, ...
 T_21_max_contact_1_p2, T_21_sum_contact_1_p2, T_22_max_contact_1_p2, T_22_sum_contact_1_p2] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p2, N, r_N, c_N_1);

                       
% 2. c = 0.4
c_N_2 = 0.4*log(log(N));

[T_11_max_contact_2_p1, T_11_sum_contact_2_p1, T_12_max_contact_2_p1, T_12_sum_contact_2_p1, ...
 T_21_max_contact_2_p1, T_21_sum_contact_2_p1, T_22_max_contact_2_p1, T_22_sum_contact_2_p1] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p1, N, r_N, c_N_2);

[T_11_max_contact_2_p2, T_11_sum_contact_2_p2, T_12_max_contact_2_p2, T_12_sum_contact_2_p2, ...
 T_21_max_contact_2_p2, T_21_sum_contact_2_p2, T_22_max_contact_2_p2, T_22_sum_contact_2_p2] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p2, N, r_N, c_N_2);
                          
% 3. c = 0.5
c_N_3 = 0.5*log(log(N));

[T_11_max_contact_3_p1, T_11_sum_contact_3_p1, T_12_max_contact_3_p1, T_12_sum_contact_3_p1, ...
 T_21_max_contact_3_p1, T_21_sum_contact_3_p1, T_22_max_contact_3_p1, T_22_sum_contact_3_p1] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p1, N, r_N, c_N_3);

[T_11_max_contact_3_p2, T_11_sum_contact_3_p2, T_12_max_contact_3_p2, T_12_sum_contact_3_p2, ...
 T_21_max_contact_3_p2, T_21_sum_contact_3_p2, T_22_max_contact_3_p2, T_22_sum_contact_3_p2] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p2, N, r_N, c_N_3);
%                                  
% 4. c = 0.6
c_N_4 = 0.6*log(log(N));

[T_11_max_contact_4_p1, T_11_sum_contact_4_p1, T_12_max_contact_4_p1, T_12_sum_contact_4_p1, ...
 T_21_max_contact_4_p1, T_21_sum_contact_4_p1, T_22_max_contact_4_p1, T_22_sum_contact_4_p1] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p1, N, r_N, c_N_4);

[T_11_max_contact_4_p2, T_11_sum_contact_4_p2, T_12_max_contact_4_p2, T_12_sum_contact_4_p2, ...
 T_21_max_contact_4_p2, T_21_sum_contact_4_p2, T_22_max_contact_4_p2, T_22_sum_contact_4_p2] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p2, N, r_N, c_N_4);
                                  
% 5. c = 0.7
c_N_5 = 0.7*log(log(N));

[T_11_max_contact_5_p1, T_11_sum_contact_5_p1, T_12_max_contact_5_p1, T_12_sum_contact_5_p1, ...
 T_21_max_contact_5_p1, T_21_sum_contact_5_p1, T_22_max_contact_5_p1, T_22_sum_contact_5_p1] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p1, N, r_N, c_N_5);


[T_11_max_contact_5_p2, T_11_sum_contact_5_p2, T_12_max_contact_5_p2, T_12_sum_contact_5_p2, ...
 T_21_max_contact_5_p2, T_21_sum_contact_5_p2, T_22_max_contact_5_p2, T_22_sum_contact_5_p2] = ...
            TGI_Estimation(D_11,D_12,D_21,D_22_collection,b_D_11_recentered,b_D_12_recentered, ...
                           b_D_21_recentered, b_D_22_collection_recentered, p2, N, r_N, c_N_5);                       
                
%--------- Numerical Delta Method ----------------------%
% We follow Hong and Li(2018)
% Refer to numerical_delta_method.m for the detail.

% About tuning Params
% We compare five epsilon parameters

%1 sqrt(log(r_N))/sqrt(r_N);
epsilon_1 = r_N^(-1/128);

% Bootstrap Statistics ----------------------%

% (1,1)
btsp_phi_dist_11_max_1_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_1,r_N,p1,'max',1);
btsp_phi_dist_11_sum_1_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_1,r_N,p1,'sum',1);

btsp_phi_dist_11_max_1_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_1,r_N,p2,'max',1);
btsp_phi_dist_11_sum_1_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_1,r_N,p2,'sum',1);

btsp_phi_dist_11_max_1_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_1,r_N,p2,'max',2);
btsp_phi_dist_11_sum_1_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_1,r_N,p2,'sum',2);

% (1,2)
btsp_phi_dist_12_max_1_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_1,r_N,p1,'max',1);
btsp_phi_dist_12_sum_1_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_1,r_N,p1,'sum',1);

btsp_phi_dist_12_max_1_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_1,r_N,p2,'max',1);
btsp_phi_dist_12_sum_1_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_1,r_N,p2,'sum',1);

btsp_phi_dist_12_max_1_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_1,r_N,p2,'max',2);
btsp_phi_dist_12_sum_1_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_1,r_N,p2,'sum',2);

% (2,1)
btsp_phi_dist_21_max_1_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_1,r_N,p1,'max',1);
btsp_phi_dist_21_sum_1_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_1,r_N,p1,'sum',1);

btsp_phi_dist_21_max_1_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_1,r_N,p2,'max',1);
btsp_phi_dist_21_sum_1_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_1,r_N,p2,'sum',1);

btsp_phi_dist_21_max_1_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_1,r_N,p2,'max',2);
btsp_phi_dist_21_sum_1_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_1,r_N,p2,'sum',2);

% (2,2)
btsp_phi_dist_22_max_1_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_1,r_N,p1,'max',1);
btsp_phi_dist_22_sum_1_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_1,r_N,p1,'sum',1);

btsp_phi_dist_22_max_1_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_1,r_N,p2,'max',1);
btsp_phi_dist_22_sum_1_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_1,r_N,p2,'sum',1);

btsp_phi_dist_22_max_1_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_1,r_N,p2,'max',2);
btsp_phi_dist_22_sum_1_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_1,r_N,p2,'sum',2);

% %2. r_N^(-1/64)
epsilon_2 = r_N^(-1/64);
% 
% % Bootstrap Statistics ----------------------%
% 
% % (1,1)
btsp_phi_dist_11_max_2_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_2,r_N,p1,'max',1);
btsp_phi_dist_11_sum_2_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_2,r_N,p1,'sum',1);

btsp_phi_dist_11_max_2_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_2,r_N,p2,'max',1);
btsp_phi_dist_11_sum_2_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_2,r_N,p2,'sum',1);

btsp_phi_dist_11_max_2_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_2,r_N,p2,'max',2);
btsp_phi_dist_11_sum_2_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_2,r_N,p2,'sum',2);

% % (1,2)
btsp_phi_dist_12_max_2_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_2,r_N,p1,'max',1);
btsp_phi_dist_12_sum_2_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_2,r_N,p1,'sum',1);

btsp_phi_dist_12_max_2_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_2,r_N,p2,'max',1);
btsp_phi_dist_12_sum_2_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_2,r_N,p2,'sum',1);

btsp_phi_dist_12_max_2_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_2,r_N,p2,'max',2);
btsp_phi_dist_12_sum_2_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_2,r_N,p2,'sum',2);

% % (2,1)
btsp_phi_dist_21_max_2_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_2,r_N,p1,'max',1);
btsp_phi_dist_21_sum_2_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_2,r_N,p1,'sum',1);

btsp_phi_dist_21_max_2_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_2,r_N,p2,'max',1);
btsp_phi_dist_21_sum_2_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_2,r_N,p2,'sum',1);

btsp_phi_dist_21_max_2_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_2,r_N,p2,'max',2);
btsp_phi_dist_21_sum_2_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_2,r_N,p2,'sum',2);

% % (2,2)
btsp_phi_dist_22_max_2_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_2,r_N,p1,'max',1);
btsp_phi_dist_22_sum_2_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_2,r_N,p1,'sum',1);

btsp_phi_dist_22_max_2_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_2,r_N,p2,'max',1);
btsp_phi_dist_22_sum_2_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_2,r_N,p2,'sum',1);

btsp_phi_dist_22_max_2_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_2,r_N,p2,'max',2);
btsp_phi_dist_22_sum_2_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_2,r_N,p2,'sum',2);

% %3. r_N^(-1/32)

epsilon_3 = r_N^(-1/32);
% 
% % Bootstrap Statistics ----------------------%
% 
% % (1,1)
btsp_phi_dist_11_max_3_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_3,r_N,p1,'max',1);
btsp_phi_dist_11_sum_3_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_3,r_N,p1,'sum',1);

btsp_phi_dist_11_max_3_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_3,r_N,p2,'max',1);
btsp_phi_dist_11_sum_3_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_3,r_N,p2,'sum',1);

btsp_phi_dist_11_max_3_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_3,r_N,p2,'max',2);
btsp_phi_dist_11_sum_3_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_3,r_N,p2,'sum',2);

% % (1,2)
btsp_phi_dist_12_max_3_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_3,r_N,p1,'max',1);
btsp_phi_dist_12_sum_3_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_3,r_N,p1,'sum',1);

btsp_phi_dist_12_max_3_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_3,r_N,p2,'max',1);
btsp_phi_dist_12_sum_3_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_3,r_N,p2,'sum',1);

btsp_phi_dist_12_max_3_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_3,r_N,p2,'max',2);
btsp_phi_dist_12_sum_3_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_3,r_N,p2,'sum',2);

% % (2,1)
btsp_phi_dist_21_max_3_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_3,r_N,p1,'max',1);
btsp_phi_dist_21_sum_3_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_3,r_N,p1,'sum',1);

btsp_phi_dist_21_max_3_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_3,r_N,p2,'max',1);
btsp_phi_dist_21_sum_3_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_3,r_N,p2,'sum',1);

btsp_phi_dist_21_max_3_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_3,r_N,p2,'max',2);
btsp_phi_dist_21_sum_3_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_3,r_N,p2,'sum',2);

% % (2,2)
btsp_phi_dist_22_max_3_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_3,r_N,p1,'max',1);
btsp_phi_dist_22_sum_3_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_3,r_N,p1,'sum',1);

btsp_phi_dist_22_max_3_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_3,r_N,p2,'max',1);
btsp_phi_dist_22_sum_3_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_3,r_N,p2,'sum',1);

btsp_phi_dist_22_max_3_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_3,r_N,p2,'max',2);
btsp_phi_dist_22_sum_3_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_3,r_N,p2,'sum',2);

% % 4. r_N^(-1/16)
% 
epsilon_4 = r_N^(-1/16);
% 
% % Bootstrap Statistics ----------------------%

% (1,1)
btsp_phi_dist_11_max_4_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_4,r_N,p1,'max',1);
btsp_phi_dist_11_sum_4_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_4,r_N,p1,'sum',1);

btsp_phi_dist_11_max_4_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_4,r_N,p2,'max',1);
btsp_phi_dist_11_sum_4_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_4,r_N,p2,'sum',1);

btsp_phi_dist_11_max_4_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_4,r_N,p2,'max',2);
btsp_phi_dist_11_sum_4_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_4,r_N,p2,'sum',2);

% (1,2)
btsp_phi_dist_12_max_4_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_4,r_N,p1,'max',1);
btsp_phi_dist_12_sum_4_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_4,r_N,p1,'sum',1);

btsp_phi_dist_12_max_4_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_4,r_N,p2,'max',1);
btsp_phi_dist_12_sum_4_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_4,r_N,p2,'sum',1);

btsp_phi_dist_12_max_4_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_4,r_N,p2,'max',2);
btsp_phi_dist_12_sum_4_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_4,r_N,p2,'sum',2);

% (2,1)
btsp_phi_dist_21_max_4_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_4,r_N,p1,'max',1);
btsp_phi_dist_21_sum_4_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_4,r_N,p1,'sum',1);

btsp_phi_dist_21_max_4_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_4,r_N,p2,'max',1);
btsp_phi_dist_21_sum_4_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_4,r_N,p2,'sum',1);

btsp_phi_dist_21_max_4_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_4,r_N,p2,'max',2);
btsp_phi_dist_21_sum_4_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_4,r_N,p2,'sum',2);

% (2,2)
btsp_phi_dist_22_max_4_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_4,r_N,p1,'max',1);
btsp_phi_dist_22_sum_4_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_4,r_N,p1,'sum',1);

btsp_phi_dist_22_max_4_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_4,r_N,p2,'max',1);
btsp_phi_dist_22_sum_4_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_4,r_N,p2,'sum',1);

btsp_phi_dist_22_max_4_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_4,r_N,p2,'max',2);
btsp_phi_dist_22_sum_4_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_4,r_N,p2,'sum',2);

% 5. r_N^(-1/8)

epsilon_5 = r_N^(-1/8);

% Bootstrap Statistics ----------------------%

% (1,1)
btsp_phi_dist_11_max_5_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_5,r_N,p1,'max',1);
btsp_phi_dist_11_sum_5_p1        = numerical_delta_method(D_11,b_D_11_recentered,epsilon_5,r_N,p1,'sum',1);

btsp_phi_dist_11_max_5_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_5,r_N,p2,'max',1);
btsp_phi_dist_11_sum_5_p2_first  = numerical_delta_method(D_11,b_D_11_recentered,epsilon_5,r_N,p2,'sum',1);

btsp_phi_dist_11_max_5_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_5,r_N,p2,'max',2);
btsp_phi_dist_11_sum_5_p2_second = numerical_delta_method(D_11,b_D_11_recentered,epsilon_5,r_N,p2,'sum',2);

% (1,2)
btsp_phi_dist_12_max_5_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_5,r_N,p1,'max',1);
btsp_phi_dist_12_sum_5_p1        = numerical_delta_method(D_12,b_D_12_recentered,epsilon_5,r_N,p1,'sum',1);

btsp_phi_dist_12_max_5_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_5,r_N,p2,'max',1);
btsp_phi_dist_12_sum_5_p2_first  = numerical_delta_method(D_12,b_D_12_recentered,epsilon_5,r_N,p2,'sum',1);

btsp_phi_dist_12_max_5_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_5,r_N,p2,'max',2);
btsp_phi_dist_12_sum_5_p2_second = numerical_delta_method(D_12,b_D_12_recentered,epsilon_5,r_N,p2,'sum',2);

% (2,1)
btsp_phi_dist_21_max_5_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_5,r_N,p1,'max',1);
btsp_phi_dist_21_sum_5_p1        = numerical_delta_method(D_21,b_D_21_recentered,epsilon_5,r_N,p1,'sum',1);

btsp_phi_dist_21_max_5_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_5,r_N,p2,'max',1);
btsp_phi_dist_21_sum_5_p2_first  = numerical_delta_method(D_21,b_D_21_recentered,epsilon_5,r_N,p2,'sum',1);

btsp_phi_dist_21_max_5_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_5,r_N,p2,'max',2);
btsp_phi_dist_21_sum_5_p2_second = numerical_delta_method(D_21,b_D_21_recentered,epsilon_5,r_N,p2,'sum',2);

% (2,2)
btsp_phi_dist_22_max_5_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_5,r_N,p1,'max',1);
btsp_phi_dist_22_sum_5_p1        = numerical_delta_method(D_22,b_D_22_recentered,epsilon_5,r_N,p1,'sum',1);

btsp_phi_dist_22_max_5_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_5,r_N,p2,'max',1);
btsp_phi_dist_22_sum_5_p2_first  = numerical_delta_method(D_22,b_D_22_recentered,epsilon_5,r_N,p2,'sum',1);

btsp_phi_dist_22_max_5_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_5,r_N,p2,'max',2);
btsp_phi_dist_22_sum_5_p2_second = numerical_delta_method(D_22,b_D_22_recentered,epsilon_5,r_N,p2,'sum',2);

% % Deciding rejection starts here
% 
% % (n,m) = 1,1

% p =1
critical_value_11_btsp_max_LFC_p1       = quantile(btsp_T_11_max_LFC_p1,1-alpha);
critical_value_11_btsp_sum_LFC_p1       = quantile(btsp_T_11_sum_LFC_p1,1-alpha);

critical_value_11_max_contact_1_p1 = max(quantile(T_11_max_contact_1_p1,1-alpha),eta);
critical_value_11_max_contact_2_p1 = max(quantile(T_11_max_contact_2_p1,1-alpha),eta);
critical_value_11_max_contact_3_p1 = max(quantile(T_11_max_contact_3_p1,1-alpha),eta);
critical_value_11_max_contact_4_p1 = max(quantile(T_11_max_contact_4_p1,1-alpha),eta);
critical_value_11_max_contact_5_p1 = max(quantile(T_11_max_contact_5_p1,1-alpha),eta);

critical_value_11_sum_contact_1_p1 = max(quantile(T_11_sum_contact_1_p1,1-alpha),eta);
critical_value_11_sum_contact_2_p1 = max(quantile(T_11_sum_contact_2_p1,1-alpha),eta);
critical_value_11_sum_contact_3_p1 = max(quantile(T_11_sum_contact_3_p1,1-alpha),eta);
critical_value_11_sum_contact_4_p1 = max(quantile(T_11_sum_contact_4_p1,1-alpha),eta);
critical_value_11_sum_contact_5_p1 = max(quantile(T_11_sum_contact_5_p1,1-alpha),eta);

critical_value_11_btsp_max_NDM_1_p1     = quantile(btsp_phi_dist_11_max_1_p1,1-alpha);
critical_value_11_btsp_sum_NDM_1_p1     = quantile(btsp_phi_dist_11_sum_1_p1,1-alpha);

critical_value_11_btsp_max_NDM_2_p1     = quantile(btsp_phi_dist_11_max_2_p1,1-alpha);
critical_value_11_btsp_sum_NDM_2_p1     = quantile(btsp_phi_dist_11_sum_2_p1,1-alpha);

critical_value_11_btsp_max_NDM_3_p1     = quantile(btsp_phi_dist_11_max_3_p1,1-alpha);
critical_value_11_btsp_sum_NDM_3_p1     = quantile(btsp_phi_dist_11_sum_3_p1,1-alpha);

critical_value_11_btsp_max_NDM_4_p1     = quantile(btsp_phi_dist_11_max_4_p1,1-alpha);
critical_value_11_btsp_sum_NDM_4_p1     = quantile(btsp_phi_dist_11_sum_4_p1,1-alpha);

critical_value_11_btsp_max_NDM_5_p1     = quantile(btsp_phi_dist_11_max_5_p1,1-alpha);
critical_value_11_btsp_sum_NDM_5_p1     = quantile(btsp_phi_dist_11_sum_5_p1,1-alpha);

% p = 2

critical_value_11_btsp_max_LFC_p2       = quantile(btsp_T_11_max_LFC_p2,1-alpha);
critical_value_11_btsp_sum_LFC_p2       = quantile(btsp_T_11_sum_LFC_p2,1-alpha);

critical_value_11_max_contact_1_p2 = max(quantile(T_11_max_contact_1_p2,1-alpha),eta);
critical_value_11_max_contact_2_p2 = max(quantile(T_11_max_contact_2_p2,1-alpha),eta);
critical_value_11_max_contact_3_p2 = max(quantile(T_11_max_contact_3_p2,1-alpha),eta);
critical_value_11_max_contact_4_p2 = max(quantile(T_11_max_contact_4_p2,1-alpha),eta);
critical_value_11_max_contact_5_p2 = max(quantile(T_11_max_contact_5_p2,1-alpha),eta);

critical_value_11_sum_contact_1_p2 = max(quantile(T_11_sum_contact_1_p2,1-alpha),eta);
critical_value_11_sum_contact_2_p2 = max(quantile(T_11_sum_contact_2_p2,1-alpha),eta);
critical_value_11_sum_contact_3_p2 = max(quantile(T_11_sum_contact_3_p2,1-alpha),eta);
critical_value_11_sum_contact_4_p2 = max(quantile(T_11_sum_contact_4_p2,1-alpha),eta);
critical_value_11_sum_contact_5_p2 = max(quantile(T_11_sum_contact_5_p2,1-alpha),eta);

% first

critical_value_11_btsp_max_NDM_1_p2_first     = quantile(btsp_phi_dist_11_max_1_p2_first,1-alpha);
critical_value_11_btsp_sum_NDM_1_p2_first     = quantile(btsp_phi_dist_11_sum_1_p2_first,1-alpha);

critical_value_11_btsp_max_NDM_2_p2_first     = quantile(btsp_phi_dist_11_max_2_p2_first,1-alpha);
critical_value_11_btsp_sum_NDM_2_p2_first     = quantile(btsp_phi_dist_11_sum_2_p2_first,1-alpha);

critical_value_11_btsp_max_NDM_3_p2_first     = quantile(btsp_phi_dist_11_max_3_p2_first,1-alpha);
critical_value_11_btsp_sum_NDM_3_p2_first     = quantile(btsp_phi_dist_11_sum_3_p2_first,1-alpha);

critical_value_11_btsp_max_NDM_4_p2_first     = quantile(btsp_phi_dist_11_max_4_p2_first,1-alpha);
critical_value_11_btsp_sum_NDM_4_p2_first     = quantile(btsp_phi_dist_11_sum_4_p2_first,1-alpha);

critical_value_11_btsp_max_NDM_5_p2_first     = quantile(btsp_phi_dist_11_max_5_p2_first,1-alpha);
critical_value_11_btsp_sum_NDM_5_p2_first     = quantile(btsp_phi_dist_11_sum_5_p2_first,1-alpha);

% second

critical_value_11_btsp_max_NDM_1_p2_second    = quantile(btsp_phi_dist_11_max_1_p2_second,1-alpha);
critical_value_11_btsp_sum_NDM_1_p2_second    = quantile(btsp_phi_dist_11_sum_1_p2_second,1-alpha);

critical_value_11_btsp_max_NDM_2_p2_second    = quantile(btsp_phi_dist_11_max_2_p2_second,1-alpha);
critical_value_11_btsp_sum_NDM_2_p2_second    = quantile(btsp_phi_dist_11_sum_2_p2_second,1-alpha);

critical_value_11_btsp_max_NDM_3_p2_second    = quantile(btsp_phi_dist_11_max_3_p2_second,1-alpha);
critical_value_11_btsp_sum_NDM_3_p2_second    = quantile(btsp_phi_dist_11_sum_3_p2_second,1-alpha);

critical_value_11_btsp_max_NDM_4_p2_second    = quantile(btsp_phi_dist_11_max_4_p2_second,1-alpha);
critical_value_11_btsp_sum_NDM_4_p2_second    = quantile(btsp_phi_dist_11_sum_4_p2_second,1-alpha);

critical_value_11_btsp_max_NDM_5_p2_second    = quantile(btsp_phi_dist_11_max_5_p2_second,1-alpha);
critical_value_11_btsp_sum_NDM_5_p2_second    = quantile(btsp_phi_dist_11_sum_5_p2_second,1-alpha);

% (n,m) = 1,2

% p =1
critical_value_12_btsp_max_LFC_p1       = quantile(btsp_T_12_max_LFC_p1,1-alpha);
critical_value_12_btsp_sum_LFC_p1       = quantile(btsp_T_12_sum_LFC_p1,1-alpha);

critical_value_12_max_contact_1_p1 = max(quantile(T_12_max_contact_1_p1,1-alpha),eta);
critical_value_12_max_contact_2_p1 = max(quantile(T_12_max_contact_2_p1,1-alpha),eta);
critical_value_12_max_contact_3_p1 = max(quantile(T_12_max_contact_3_p1,1-alpha),eta);
critical_value_12_max_contact_4_p1 = max(quantile(T_12_max_contact_4_p1,1-alpha),eta);
critical_value_12_max_contact_5_p1 = max(quantile(T_12_max_contact_5_p1,1-alpha),eta);

critical_value_12_sum_contact_1_p1 = max(quantile(T_12_sum_contact_1_p1,1-alpha),eta);
critical_value_12_sum_contact_2_p1 = max(quantile(T_12_sum_contact_2_p1,1-alpha),eta);
critical_value_12_sum_contact_3_p1 = max(quantile(T_12_sum_contact_3_p1,1-alpha),eta);
critical_value_12_sum_contact_4_p1 = max(quantile(T_12_sum_contact_4_p1,1-alpha),eta);
critical_value_12_sum_contact_5_p1 = max(quantile(T_12_sum_contact_5_p1,1-alpha),eta);

critical_value_12_btsp_max_NDM_1_p1     = quantile(btsp_phi_dist_12_max_1_p1,1-alpha);
critical_value_12_btsp_sum_NDM_1_p1     = quantile(btsp_phi_dist_12_sum_1_p1,1-alpha);
 
critical_value_12_btsp_max_NDM_2_p1     = quantile(btsp_phi_dist_12_max_2_p1,1-alpha);
critical_value_12_btsp_sum_NDM_2_p1     = quantile(btsp_phi_dist_12_sum_2_p1,1-alpha);

critical_value_12_btsp_max_NDM_3_p1     = quantile(btsp_phi_dist_12_max_3_p1,1-alpha);
critical_value_12_btsp_sum_NDM_3_p1     = quantile(btsp_phi_dist_12_sum_3_p1,1-alpha);

critical_value_12_btsp_max_NDM_4_p1     = quantile(btsp_phi_dist_12_max_4_p1,1-alpha);
critical_value_12_btsp_sum_NDM_4_p1     = quantile(btsp_phi_dist_12_sum_4_p1,1-alpha);

critical_value_12_btsp_max_NDM_5_p1     = quantile(btsp_phi_dist_12_max_5_p1,1-alpha);
critical_value_12_btsp_sum_NDM_5_p1     = quantile(btsp_phi_dist_12_sum_5_p1,1-alpha);

% p = 2

critical_value_12_btsp_max_LFC_p2       = quantile(btsp_T_12_max_LFC_p2,1-alpha);
critical_value_12_btsp_sum_LFC_p2       = quantile(btsp_T_12_sum_LFC_p2,1-alpha);

critical_value_12_max_contact_1_p2 = max(quantile(T_12_max_contact_1_p2,1-alpha),eta);
critical_value_12_max_contact_2_p2 = max(quantile(T_12_max_contact_2_p2,1-alpha),eta);
critical_value_12_max_contact_3_p2 = max(quantile(T_12_max_contact_3_p2,1-alpha),eta);
critical_value_12_max_contact_4_p2 = max(quantile(T_12_max_contact_4_p2,1-alpha),eta);
critical_value_12_max_contact_5_p2 = max(quantile(T_12_max_contact_5_p2,1-alpha),eta);

critical_value_12_sum_contact_1_p2 = max(quantile(T_12_sum_contact_1_p2,1-alpha),eta);
critical_value_12_sum_contact_2_p2 = max(quantile(T_12_sum_contact_2_p2,1-alpha),eta);
critical_value_12_sum_contact_3_p2 = max(quantile(T_12_sum_contact_3_p2,1-alpha),eta);
critical_value_12_sum_contact_4_p2 = max(quantile(T_12_sum_contact_4_p2,1-alpha),eta);
critical_value_12_sum_contact_5_p2 = max(quantile(T_12_sum_contact_5_p2,1-alpha),eta);

% first

critical_value_12_btsp_max_NDM_1_p2_first     = quantile(btsp_phi_dist_12_max_1_p2_first,1-alpha);
critical_value_12_btsp_sum_NDM_1_p2_first     = quantile(btsp_phi_dist_12_sum_1_p2_first,1-alpha);

critical_value_12_btsp_max_NDM_2_p2_first     = quantile(btsp_phi_dist_12_max_2_p2_first,1-alpha);
critical_value_12_btsp_sum_NDM_2_p2_first     = quantile(btsp_phi_dist_12_sum_2_p2_first,1-alpha);

critical_value_12_btsp_max_NDM_3_p2_first     = quantile(btsp_phi_dist_12_max_3_p2_first,1-alpha);
critical_value_12_btsp_sum_NDM_3_p2_first     = quantile(btsp_phi_dist_12_sum_3_p2_first,1-alpha);

critical_value_12_btsp_max_NDM_4_p2_first     = quantile(btsp_phi_dist_12_max_4_p2_first,1-alpha);
critical_value_12_btsp_sum_NDM_4_p2_first     = quantile(btsp_phi_dist_12_sum_4_p2_first,1-alpha);

critical_value_12_btsp_max_NDM_5_p2_first     = quantile(btsp_phi_dist_12_max_5_p2_first,1-alpha);
critical_value_12_btsp_sum_NDM_5_p2_first     = quantile(btsp_phi_dist_12_sum_5_p2_first,1-alpha);

% second

critical_value_12_btsp_max_NDM_1_p2_second    = quantile(btsp_phi_dist_12_max_1_p2_second,1-alpha);
critical_value_12_btsp_sum_NDM_1_p2_second    = quantile(btsp_phi_dist_12_sum_1_p2_second,1-alpha);
% 
critical_value_12_btsp_max_NDM_2_p2_second    = quantile(btsp_phi_dist_12_max_2_p2_second,1-alpha);
critical_value_12_btsp_sum_NDM_2_p2_second    = quantile(btsp_phi_dist_12_sum_2_p2_second,1-alpha);

critical_value_12_btsp_max_NDM_3_p2_second    = quantile(btsp_phi_dist_12_max_3_p2_second,1-alpha);
critical_value_12_btsp_sum_NDM_3_p2_second    = quantile(btsp_phi_dist_12_sum_3_p2_second,1-alpha);

critical_value_12_btsp_max_NDM_4_p2_second    = quantile(btsp_phi_dist_12_max_4_p2_second,1-alpha);
critical_value_12_btsp_sum_NDM_4_p2_second    = quantile(btsp_phi_dist_12_sum_4_p2_second,1-alpha);

critical_value_12_btsp_max_NDM_5_p2_second    = quantile(btsp_phi_dist_12_max_5_p2_second,1-alpha);
critical_value_12_btsp_sum_NDM_5_p2_second    = quantile(btsp_phi_dist_12_sum_5_p2_second,1-alpha);

% (n,m) = 2,1

% p =1
critical_value_21_btsp_max_LFC_p1       = quantile(btsp_T_21_max_LFC_p1,1-alpha);
critical_value_21_btsp_sum_LFC_p1       = quantile(btsp_T_21_sum_LFC_p1,1-alpha);

critical_value_21_max_contact_1_p1 = max(quantile(T_21_max_contact_1_p1,1-alpha),eta);
critical_value_21_max_contact_2_p1 = max(quantile(T_21_max_contact_2_p1,1-alpha),eta);
critical_value_21_max_contact_3_p1 = max(quantile(T_21_max_contact_3_p1,1-alpha),eta);
critical_value_21_max_contact_4_p1 = max(quantile(T_21_max_contact_4_p1,1-alpha),eta);
critical_value_21_max_contact_5_p1 = max(quantile(T_21_max_contact_5_p1,1-alpha),eta);

critical_value_21_sum_contact_1_p1 = max(quantile(T_21_sum_contact_1_p1,1-alpha),eta);
critical_value_21_sum_contact_2_p1 = max(quantile(T_21_sum_contact_2_p1,1-alpha),eta);
critical_value_21_sum_contact_3_p1 = max(quantile(T_21_sum_contact_3_p1,1-alpha),eta);
critical_value_21_sum_contact_4_p1 = max(quantile(T_21_sum_contact_4_p1,1-alpha),eta);
critical_value_21_sum_contact_5_p1 = max(quantile(T_21_sum_contact_5_p1,1-alpha),eta);

critical_value_21_btsp_max_NDM_1_p1     = quantile(btsp_phi_dist_21_max_1_p1,1-alpha);
critical_value_21_btsp_sum_NDM_1_p1     = quantile(btsp_phi_dist_21_sum_1_p1,1-alpha);

critical_value_21_btsp_max_NDM_2_p1     = quantile(btsp_phi_dist_21_max_2_p1,1-alpha);
critical_value_21_btsp_sum_NDM_2_p1     = quantile(btsp_phi_dist_21_sum_2_p1,1-alpha);

critical_value_21_btsp_max_NDM_3_p1     = quantile(btsp_phi_dist_21_max_3_p1,1-alpha);
critical_value_21_btsp_sum_NDM_3_p1     = quantile(btsp_phi_dist_21_sum_3_p1,1-alpha);

critical_value_21_btsp_max_NDM_4_p1     = quantile(btsp_phi_dist_21_max_4_p1,1-alpha);
critical_value_21_btsp_sum_NDM_4_p1     = quantile(btsp_phi_dist_21_sum_4_p1,1-alpha);

critical_value_21_btsp_max_NDM_5_p1     = quantile(btsp_phi_dist_21_max_5_p1,1-alpha);
critical_value_21_btsp_sum_NDM_5_p1     = quantile(btsp_phi_dist_21_sum_5_p1,1-alpha);

% p = 2

critical_value_21_btsp_max_LFC_p2       = quantile(btsp_T_21_max_LFC_p2,1-alpha);
critical_value_21_btsp_sum_LFC_p2       = quantile(btsp_T_21_sum_LFC_p2,1-alpha);

critical_value_21_max_contact_1_p2 = max(quantile(T_21_max_contact_1_p2,1-alpha),eta);
critical_value_21_max_contact_2_p2 = max(quantile(T_21_max_contact_2_p2,1-alpha),eta);
critical_value_21_max_contact_3_p2 = max(quantile(T_21_max_contact_3_p2,1-alpha),eta);
critical_value_21_max_contact_4_p2 = max(quantile(T_21_max_contact_4_p2,1-alpha),eta);
critical_value_21_max_contact_5_p2 = max(quantile(T_21_max_contact_5_p2,1-alpha),eta);

critical_value_21_sum_contact_1_p2 = max(quantile(T_21_sum_contact_1_p2,1-alpha),eta);
critical_value_21_sum_contact_2_p2 = max(quantile(T_21_sum_contact_2_p2,1-alpha),eta);
critical_value_21_sum_contact_3_p2 = max(quantile(T_21_sum_contact_3_p2,1-alpha),eta);
critical_value_21_sum_contact_4_p2 = max(quantile(T_21_sum_contact_4_p2,1-alpha),eta);
critical_value_21_sum_contact_5_p2 = max(quantile(T_21_sum_contact_5_p2,1-alpha),eta);

% first

critical_value_21_btsp_max_NDM_1_p2_first     = quantile(btsp_phi_dist_21_max_1_p2_first,1-alpha);
critical_value_21_btsp_sum_NDM_1_p2_first     = quantile(btsp_phi_dist_21_sum_1_p2_first,1-alpha);

critical_value_21_btsp_max_NDM_2_p2_first     = quantile(btsp_phi_dist_21_max_2_p2_first,1-alpha);
critical_value_21_btsp_sum_NDM_2_p2_first     = quantile(btsp_phi_dist_21_sum_2_p2_first,1-alpha);

critical_value_21_btsp_max_NDM_3_p2_first     = quantile(btsp_phi_dist_21_max_3_p2_first,1-alpha);
critical_value_21_btsp_sum_NDM_3_p2_first     = quantile(btsp_phi_dist_21_sum_3_p2_first,1-alpha);

critical_value_21_btsp_max_NDM_4_p2_first     = quantile(btsp_phi_dist_21_max_4_p2_first,1-alpha);
critical_value_21_btsp_sum_NDM_4_p2_first     = quantile(btsp_phi_dist_21_sum_4_p2_first,1-alpha);

critical_value_21_btsp_max_NDM_5_p2_first     = quantile(btsp_phi_dist_21_max_5_p2_first,1-alpha);
critical_value_21_btsp_sum_NDM_5_p2_first     = quantile(btsp_phi_dist_21_sum_5_p2_first,1-alpha);

% second

critical_value_21_btsp_max_NDM_1_p2_second    = quantile(btsp_phi_dist_21_max_1_p2_second,1-alpha);
critical_value_21_btsp_sum_NDM_1_p2_second    = quantile(btsp_phi_dist_21_sum_1_p2_second,1-alpha);

critical_value_21_btsp_max_NDM_2_p2_second    = quantile(btsp_phi_dist_21_max_2_p2_second,1-alpha);
critical_value_21_btsp_sum_NDM_2_p2_second    = quantile(btsp_phi_dist_21_sum_2_p2_second,1-alpha);

critical_value_21_btsp_max_NDM_3_p2_second    = quantile(btsp_phi_dist_21_max_3_p2_second,1-alpha);
critical_value_21_btsp_sum_NDM_3_p2_second    = quantile(btsp_phi_dist_21_sum_3_p2_second,1-alpha);

critical_value_21_btsp_max_NDM_4_p2_second    = quantile(btsp_phi_dist_21_max_4_p2_second,1-alpha);
critical_value_21_btsp_sum_NDM_4_p2_second    = quantile(btsp_phi_dist_21_sum_4_p2_second,1-alpha);

critical_value_21_btsp_max_NDM_5_p2_second    = quantile(btsp_phi_dist_21_max_5_p2_second,1-alpha);
critical_value_21_btsp_sum_NDM_5_p2_second    = quantile(btsp_phi_dist_21_sum_5_p2_second,1-alpha);

% (n,m) = 2,2
% p = 1

% p =1
critical_value_22_btsp_max_LFC_p1       = quantile(btsp_T_22_max_LFC_p1,1-alpha);
critical_value_22_btsp_sum_LFC_p1       = quantile(btsp_T_22_sum_LFC_p1,1-alpha);

critical_value_22_max_contact_1_p1 = max(quantile(T_22_max_contact_1_p1,1-alpha),eta);
critical_value_22_max_contact_2_p1 = max(quantile(T_22_max_contact_2_p1,1-alpha),eta);
critical_value_22_max_contact_3_p1 = max(quantile(T_22_max_contact_3_p1,1-alpha),eta);
critical_value_22_max_contact_4_p1 = max(quantile(T_22_max_contact_4_p1,1-alpha),eta);
critical_value_22_max_contact_5_p1 = max(quantile(T_22_max_contact_5_p1,1-alpha),eta);

critical_value_22_sum_contact_1_p1 = max(quantile(T_22_sum_contact_1_p1,1-alpha),eta);
critical_value_22_sum_contact_2_p1 = max(quantile(T_22_sum_contact_2_p1,1-alpha),eta);
critical_value_22_sum_contact_3_p1 = max(quantile(T_22_sum_contact_3_p1,1-alpha),eta);
critical_value_22_sum_contact_4_p1 = max(quantile(T_22_sum_contact_4_p1,1-alpha),eta);
critical_value_22_sum_contact_5_p1 = max(quantile(T_22_sum_contact_5_p1,1-alpha),eta);

critical_value_22_btsp_max_NDM_1_p1     = quantile(btsp_phi_dist_22_max_1_p1,1-alpha);
critical_value_22_btsp_sum_NDM_1_p1     = quantile(btsp_phi_dist_22_sum_1_p1,1-alpha);
% 
critical_value_22_btsp_max_NDM_2_p1     = quantile(btsp_phi_dist_22_max_2_p1,1-alpha);
critical_value_22_btsp_sum_NDM_2_p1     = quantile(btsp_phi_dist_22_sum_2_p1,1-alpha);

critical_value_22_btsp_max_NDM_3_p1     = quantile(btsp_phi_dist_22_max_3_p1,1-alpha);
critical_value_22_btsp_sum_NDM_3_p1     = quantile(btsp_phi_dist_22_sum_3_p1,1-alpha);

critical_value_22_btsp_max_NDM_4_p1     = quantile(btsp_phi_dist_22_max_4_p1,1-alpha);
critical_value_22_btsp_sum_NDM_4_p1     = quantile(btsp_phi_dist_22_sum_4_p1,1-alpha);

critical_value_22_btsp_max_NDM_5_p1     = quantile(btsp_phi_dist_22_max_5_p1,1-alpha);
critical_value_22_btsp_sum_NDM_5_p1     = quantile(btsp_phi_dist_22_sum_5_p1,1-alpha);

% p = 2

critical_value_22_btsp_max_LFC_p2       = quantile(btsp_T_22_max_LFC_p2,1-alpha);
critical_value_22_btsp_sum_LFC_p2       = quantile(btsp_T_22_sum_LFC_p2,1-alpha);

critical_value_22_max_contact_1_p2 = max(quantile(T_22_max_contact_1_p2,1-alpha),eta);
critical_value_22_max_contact_2_p2 = max(quantile(T_22_max_contact_2_p2,1-alpha),eta);
critical_value_22_max_contact_3_p2 = max(quantile(T_22_max_contact_3_p2,1-alpha),eta);
critical_value_22_max_contact_4_p2 = max(quantile(T_22_max_contact_4_p2,1-alpha),eta);
critical_value_22_max_contact_5_p2 = max(quantile(T_22_max_contact_5_p2,1-alpha),eta);

critical_value_22_sum_contact_1_p2 = max(quantile(T_22_sum_contact_1_p2,1-alpha),eta);
critical_value_22_sum_contact_2_p2 = max(quantile(T_22_sum_contact_2_p2,1-alpha),eta);
critical_value_22_sum_contact_3_p2 = max(quantile(T_22_sum_contact_3_p2,1-alpha),eta);
critical_value_22_sum_contact_4_p2 = max(quantile(T_22_sum_contact_4_p2,1-alpha),eta);
critical_value_22_sum_contact_5_p2 = max(quantile(T_22_sum_contact_5_p2,1-alpha),eta);

% first

critical_value_22_btsp_max_NDM_1_p2_first     = quantile(btsp_phi_dist_22_max_1_p2_first,1-alpha);
critical_value_22_btsp_sum_NDM_1_p2_first     = quantile(btsp_phi_dist_22_sum_1_p2_first,1-alpha);

critical_value_22_btsp_max_NDM_2_p2_first     = quantile(btsp_phi_dist_22_max_2_p2_first,1-alpha);
critical_value_22_btsp_sum_NDM_2_p2_first     = quantile(btsp_phi_dist_22_sum_2_p2_first,1-alpha);

critical_value_22_btsp_max_NDM_3_p2_first     = quantile(btsp_phi_dist_22_max_3_p2_first,1-alpha);
critical_value_22_btsp_sum_NDM_3_p2_first     = quantile(btsp_phi_dist_22_sum_3_p2_first,1-alpha);

critical_value_22_btsp_max_NDM_4_p2_first     = quantile(btsp_phi_dist_22_max_4_p2_first,1-alpha);
critical_value_22_btsp_sum_NDM_4_p2_first     = quantile(btsp_phi_dist_22_sum_4_p2_first,1-alpha);

critical_value_22_btsp_max_NDM_5_p2_first     = quantile(btsp_phi_dist_22_max_5_p2_first,1-alpha);
critical_value_22_btsp_sum_NDM_5_p2_first     = quantile(btsp_phi_dist_22_sum_5_p2_first,1-alpha);

% second

critical_value_22_btsp_max_NDM_1_p2_second    = quantile(btsp_phi_dist_22_max_1_p2_second,1-alpha);
critical_value_22_btsp_sum_NDM_1_p2_second    = quantile(btsp_phi_dist_22_sum_1_p2_second,1-alpha);

critical_value_22_btsp_max_NDM_2_p2_second    = quantile(btsp_phi_dist_22_max_2_p2_second,1-alpha);
critical_value_22_btsp_sum_NDM_2_p2_second    = quantile(btsp_phi_dist_22_sum_2_p2_second,1-alpha);

critical_value_22_btsp_max_NDM_3_p2_second    = quantile(btsp_phi_dist_22_max_3_p2_second,1-alpha);
critical_value_22_btsp_sum_NDM_3_p2_second    = quantile(btsp_phi_dist_22_sum_3_p2_second,1-alpha);

critical_value_22_btsp_max_NDM_4_p2_second    = quantile(btsp_phi_dist_22_max_4_p2_second,1-alpha);
critical_value_22_btsp_sum_NDM_4_p2_second    = quantile(btsp_phi_dist_22_sum_4_p2_second,1-alpha);

critical_value_22_btsp_max_NDM_5_p2_second    = quantile(btsp_phi_dist_22_max_5_p2_second,1-alpha);
critical_value_22_btsp_sum_NDM_5_p2_second    = quantile(btsp_phi_dist_22_sum_5_p2_second,1-alpha);

% Checking Rejection 
% % (1,1)

% p = 1

rejection_11_btsp_max_LFC_p1     = [rejection_11_btsp_max_LFC_p1, T_11_max_p1 > critical_value_11_btsp_max_LFC_p1];
rejection_11_btsp_sum_LFC_p1     = [rejection_11_btsp_sum_LFC_p1, T_11_sum_p1 > critical_value_11_btsp_sum_LFC_p1];

rejection_11_max_contact_1_p1 =  [rejection_11_max_contact_1_p1, T_11_max_p1 > critical_value_11_max_contact_1_p1];
rejection_11_sum_contact_1_p1 =  [rejection_11_sum_contact_1_p1, T_11_sum_p1 > critical_value_11_sum_contact_1_p1];

rejection_11_max_contact_2_p1 =  [rejection_11_max_contact_2_p1, T_11_max_p1 > critical_value_11_max_contact_2_p1];
rejection_11_sum_contact_2_p1 =  [rejection_11_sum_contact_2_p1, T_11_sum_p1 > critical_value_11_sum_contact_2_p1];

rejection_11_max_contact_3_p1 =  [rejection_11_max_contact_3_p1, T_11_max_p1 > critical_value_11_max_contact_3_p1];
rejection_11_sum_contact_3_p1 =  [rejection_11_sum_contact_3_p1, T_11_sum_p1 > critical_value_11_sum_contact_3_p1];

rejection_11_max_contact_4_p1 =  [rejection_11_max_contact_4_p1, T_11_max_p1 > critical_value_11_max_contact_4_p1];
rejection_11_sum_contact_4_p1 =  [rejection_11_sum_contact_4_p1, T_11_sum_p1 > critical_value_11_sum_contact_4_p1];

rejection_11_max_contact_5_p1 =  [rejection_11_max_contact_5_p1, T_11_max_p1 > critical_value_11_max_contact_5_p1];
rejection_11_sum_contact_5_p1 =  [rejection_11_sum_contact_5_p1, T_11_sum_p1 > critical_value_11_sum_contact_5_p1];

rejection_11_btsp_max_NDM_1_p1     = [rejection_11_btsp_max_NDM_1_p1, T_11_max_p1 > critical_value_11_btsp_max_NDM_1_p1];
rejection_11_btsp_sum_NDM_1_p1     = [rejection_11_btsp_sum_NDM_1_p1, T_11_sum_p1 > critical_value_11_btsp_sum_NDM_1_p1];

rejection_11_btsp_max_NDM_2_p1     = [rejection_11_btsp_max_NDM_2_p1, T_11_max_p1 > critical_value_11_btsp_max_NDM_2_p1];
rejection_11_btsp_sum_NDM_2_p1     = [rejection_11_btsp_sum_NDM_2_p1, T_11_sum_p1 > critical_value_11_btsp_sum_NDM_2_p1];

rejection_11_btsp_max_NDM_3_p1     = [rejection_11_btsp_max_NDM_3_p1, T_11_max_p1 > critical_value_11_btsp_max_NDM_3_p1];
rejection_11_btsp_sum_NDM_3_p1     = [rejection_11_btsp_sum_NDM_3_p1, T_11_sum_p1 > critical_value_11_btsp_sum_NDM_3_p1];

rejection_11_btsp_max_NDM_4_p1     = [rejection_11_btsp_max_NDM_4_p1, T_11_max_p1 > critical_value_11_btsp_max_NDM_4_p1];
rejection_11_btsp_sum_NDM_4_p1     = [rejection_11_btsp_sum_NDM_4_p1, T_11_sum_p1 > critical_value_11_btsp_sum_NDM_4_p1];

rejection_11_btsp_max_NDM_5_p1     = [rejection_11_btsp_max_NDM_5_p1, T_11_max_p1 > critical_value_11_btsp_max_NDM_5_p1];
rejection_11_btsp_sum_NDM_5_p1     = [rejection_11_btsp_sum_NDM_5_p1, T_11_sum_p1 > critical_value_11_btsp_sum_NDM_5_p1];

% p = 2

rejection_11_btsp_max_LFC_p2     = [rejection_11_btsp_max_LFC_p2, T_11_max_p2 > critical_value_11_btsp_max_LFC_p2];
rejection_11_btsp_sum_LFC_p2     = [rejection_11_btsp_sum_LFC_p2, T_11_sum_p2 > critical_value_11_btsp_sum_LFC_p2];


rejection_11_max_contact_1_p2    = [rejection_11_max_contact_1_p2, T_11_max_p2 > critical_value_11_max_contact_1_p2];
rejection_11_sum_contact_1_p2    = [rejection_11_sum_contact_1_p2, T_11_sum_p2 > critical_value_11_sum_contact_1_p2];

rejection_11_max_contact_2_p2    = [rejection_11_max_contact_2_p2, T_11_max_p2 > critical_value_11_max_contact_2_p2];
rejection_11_sum_contact_2_p2    = [rejection_11_sum_contact_2_p2, T_11_sum_p2 > critical_value_11_sum_contact_2_p2];

rejection_11_max_contact_3_p2    = [rejection_11_max_contact_3_p2, T_11_max_p2 > critical_value_11_max_contact_3_p2];
rejection_11_sum_contact_3_p2    = [rejection_11_sum_contact_3_p2, T_11_sum_p2 > critical_value_11_sum_contact_3_p2];

rejection_11_max_contact_4_p2    = [rejection_11_max_contact_4_p2, T_11_max_p2 > critical_value_11_max_contact_4_p2];
rejection_11_sum_contact_4_p2    = [rejection_11_sum_contact_4_p2, T_11_sum_p2 > critical_value_11_sum_contact_4_p2];

rejection_11_max_contact_5_p2    = [rejection_11_max_contact_5_p2, T_11_max_p2 > critical_value_11_max_contact_5_p2];
rejection_11_sum_contact_5_p2    = [rejection_11_sum_contact_5_p2, T_11_sum_p2 > critical_value_11_sum_contact_5_p2];

rejection_11_btsp_max_NDM_1_p2_first     = [rejection_11_btsp_max_NDM_1_p2_first, T_11_max_p2 > critical_value_11_btsp_max_NDM_1_p2_first];
rejection_11_btsp_sum_NDM_1_p2_first     = [rejection_11_btsp_sum_NDM_1_p2_first, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_1_p2_first];

rejection_11_btsp_max_NDM_2_p2_first     = [rejection_11_btsp_max_NDM_2_p2_first, T_11_max_p2 > critical_value_11_btsp_max_NDM_2_p2_first];
rejection_11_btsp_sum_NDM_2_p2_first     = [rejection_11_btsp_sum_NDM_2_p2_first, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_2_p2_first];

rejection_11_btsp_max_NDM_3_p2_first     = [rejection_11_btsp_max_NDM_3_p2_first, T_11_max_p2 > critical_value_11_btsp_max_NDM_3_p2_first];
rejection_11_btsp_sum_NDM_3_p2_first     = [rejection_11_btsp_sum_NDM_3_p2_first, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_3_p2_first];

rejection_11_btsp_max_NDM_4_p2_first     = [rejection_11_btsp_max_NDM_4_p2_first, T_11_max_p2 > critical_value_11_btsp_max_NDM_4_p2_first];
rejection_11_btsp_sum_NDM_4_p2_first     = [rejection_11_btsp_sum_NDM_4_p2_first, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_4_p2_first];

rejection_11_btsp_max_NDM_5_p2_first     = [rejection_11_btsp_max_NDM_5_p2_first, T_11_max_p2 > critical_value_11_btsp_max_NDM_5_p2_first];
rejection_11_btsp_sum_NDM_5_p2_first     = [rejection_11_btsp_sum_NDM_5_p2_first, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_5_p2_first];

rejection_11_btsp_max_NDM_1_p2_second    = [rejection_11_btsp_max_NDM_1_p2_second, T_11_max_p2 > critical_value_11_btsp_max_NDM_1_p2_second];
rejection_11_btsp_sum_NDM_1_p2_second    = [rejection_11_btsp_sum_NDM_1_p2_second, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_1_p2_second];

rejection_11_btsp_max_NDM_2_p2_second    = [rejection_11_btsp_max_NDM_2_p2_second, T_11_max_p2 > critical_value_11_btsp_max_NDM_2_p2_second];
rejection_11_btsp_sum_NDM_2_p2_second    = [rejection_11_btsp_sum_NDM_2_p2_second, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_2_p2_second];

rejection_11_btsp_max_NDM_3_p2_second    = [rejection_11_btsp_max_NDM_3_p2_second, T_11_max_p2 > critical_value_11_btsp_max_NDM_3_p2_second];
rejection_11_btsp_sum_NDM_3_p2_second    = [rejection_11_btsp_sum_NDM_3_p2_second, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_3_p2_second];

rejection_11_btsp_max_NDM_4_p2_second    = [rejection_11_btsp_max_NDM_4_p2_second, T_11_max_p2 > critical_value_11_btsp_max_NDM_4_p2_second];
rejection_11_btsp_sum_NDM_4_p2_second    = [rejection_11_btsp_sum_NDM_4_p2_second, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_4_p2_second];

rejection_11_btsp_max_NDM_5_p2_second    = [rejection_11_btsp_max_NDM_5_p2_second, T_11_max_p2 > critical_value_11_btsp_max_NDM_5_p2_second];
rejection_11_btsp_sum_NDM_5_p2_second    = [rejection_11_btsp_sum_NDM_5_p2_second, T_11_sum_p2 > critical_value_11_btsp_sum_NDM_5_p2_second];

% (1,2)

rejection_12_btsp_max_LFC_p1     = [rejection_12_btsp_max_LFC_p1, T_12_max_p1 > critical_value_12_btsp_max_LFC_p1];
rejection_12_btsp_sum_LFC_p1     = [rejection_12_btsp_sum_LFC_p1, T_12_sum_p1 > critical_value_12_btsp_sum_LFC_p1];

rejection_12_max_contact_1_p1 =  [rejection_12_max_contact_1_p1, T_12_max_p1 > critical_value_12_max_contact_1_p1];
rejection_12_sum_contact_1_p1 =  [rejection_12_sum_contact_1_p1, T_12_sum_p1 > critical_value_12_sum_contact_1_p1];

rejection_12_max_contact_2_p1 =  [rejection_12_max_contact_2_p1, T_12_max_p1 > critical_value_12_max_contact_2_p1];
rejection_12_sum_contact_2_p1 =  [rejection_12_sum_contact_2_p1, T_12_sum_p1 > critical_value_12_sum_contact_2_p1];

rejection_12_max_contact_3_p1 =  [rejection_12_max_contact_3_p1, T_12_max_p1 > critical_value_12_max_contact_3_p1];
rejection_12_sum_contact_3_p1 =  [rejection_12_sum_contact_3_p1, T_12_sum_p1 > critical_value_12_sum_contact_3_p1];

rejection_12_max_contact_4_p1 =  [rejection_12_max_contact_4_p1, T_12_max_p1 > critical_value_12_max_contact_4_p1];
rejection_12_sum_contact_4_p1 =  [rejection_12_sum_contact_4_p1, T_12_sum_p1 > critical_value_12_sum_contact_4_p1];

rejection_12_max_contact_5_p1 =  [rejection_12_max_contact_5_p1, T_12_max_p1 > critical_value_12_max_contact_5_p1];
rejection_12_sum_contact_5_p1 =  [rejection_12_sum_contact_5_p1, T_12_sum_p1 > critical_value_12_sum_contact_5_p1];

rejection_12_btsp_max_NDM_1_p1     = [rejection_12_btsp_max_NDM_1_p1, T_12_max_p1 > critical_value_12_btsp_max_NDM_1_p1];
rejection_12_btsp_sum_NDM_1_p1     = [rejection_12_btsp_sum_NDM_1_p1, T_12_sum_p1 > critical_value_12_btsp_sum_NDM_1_p1];

rejection_12_btsp_max_NDM_2_p1     = [rejection_12_btsp_max_NDM_2_p1, T_12_max_p1 > critical_value_12_btsp_max_NDM_2_p1];
rejection_12_btsp_sum_NDM_2_p1     = [rejection_12_btsp_sum_NDM_2_p1, T_12_sum_p1 > critical_value_12_btsp_sum_NDM_2_p1];

rejection_12_btsp_max_NDM_3_p1     = [rejection_12_btsp_max_NDM_3_p1, T_12_max_p1 > critical_value_12_btsp_max_NDM_3_p1];
rejection_12_btsp_sum_NDM_3_p1     = [rejection_12_btsp_sum_NDM_3_p1, T_12_sum_p1 > critical_value_12_btsp_sum_NDM_3_p1];

rejection_12_btsp_max_NDM_4_p1     = [rejection_12_btsp_max_NDM_4_p1, T_12_max_p1 > critical_value_12_btsp_max_NDM_4_p1];
rejection_12_btsp_sum_NDM_4_p1     = [rejection_12_btsp_sum_NDM_4_p1, T_12_sum_p1 > critical_value_12_btsp_sum_NDM_4_p1];

rejection_12_btsp_max_NDM_5_p1     = [rejection_12_btsp_max_NDM_5_p1, T_12_max_p1 > critical_value_12_btsp_max_NDM_5_p1];
rejection_12_btsp_sum_NDM_5_p1     = [rejection_12_btsp_sum_NDM_5_p1, T_12_sum_p1 > critical_value_12_btsp_sum_NDM_5_p1];

% p = 2

rejection_12_btsp_max_LFC_p2     = [rejection_12_btsp_max_LFC_p2, T_12_max_p2 > critical_value_12_btsp_max_LFC_p2];
rejection_12_btsp_sum_LFC_p2     = [rejection_12_btsp_sum_LFC_p2, T_12_sum_p2 > critical_value_12_btsp_sum_LFC_p2];

rejection_12_max_contact_1_p2    = [rejection_12_max_contact_1_p2, T_12_max_p2 > critical_value_12_max_contact_1_p2];
rejection_12_sum_contact_1_p2    = [rejection_12_sum_contact_1_p2, T_12_sum_p2 > critical_value_12_sum_contact_1_p2];

rejection_12_max_contact_2_p2    = [rejection_12_max_contact_2_p2, T_12_max_p2 > critical_value_12_max_contact_2_p2];
rejection_12_sum_contact_2_p2    = [rejection_12_sum_contact_2_p2, T_12_sum_p2 > critical_value_12_sum_contact_2_p2];

rejection_12_max_contact_3_p2    = [rejection_12_max_contact_3_p2, T_12_max_p2 > critical_value_12_max_contact_3_p2];
rejection_12_sum_contact_3_p2    = [rejection_12_sum_contact_3_p2, T_12_sum_p2 > critical_value_12_sum_contact_3_p2];

rejection_12_max_contact_4_p2    = [rejection_12_max_contact_4_p2, T_12_max_p2 > critical_value_12_max_contact_4_p2];
rejection_12_sum_contact_4_p2    = [rejection_12_sum_contact_4_p2, T_12_sum_p2 > critical_value_12_sum_contact_4_p2];

rejection_12_max_contact_5_p2    = [rejection_12_max_contact_5_p2, T_12_max_p2 > critical_value_12_max_contact_5_p2];
rejection_12_sum_contact_5_p2    = [rejection_12_sum_contact_5_p2, T_12_sum_p2 > critical_value_12_sum_contact_5_p2];

rejection_12_btsp_max_NDM_1_p2_first     = [rejection_12_btsp_max_NDM_1_p2_first, T_12_max_p2 > critical_value_12_btsp_max_NDM_1_p2_first];
rejection_12_btsp_sum_NDM_1_p2_first     = [rejection_12_btsp_sum_NDM_1_p2_first, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_1_p2_first];

rejection_12_btsp_max_NDM_2_p2_first     = [rejection_12_btsp_max_NDM_2_p2_first, T_12_max_p2 > critical_value_12_btsp_max_NDM_2_p2_first];
rejection_12_btsp_sum_NDM_2_p2_first     = [rejection_12_btsp_sum_NDM_2_p2_first, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_2_p2_first];

rejection_12_btsp_max_NDM_3_p2_first     = [rejection_12_btsp_max_NDM_3_p2_first, T_12_max_p2 > critical_value_12_btsp_max_NDM_3_p2_first];
rejection_12_btsp_sum_NDM_3_p2_first     = [rejection_12_btsp_sum_NDM_3_p2_first, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_3_p2_first];

rejection_12_btsp_max_NDM_4_p2_first     = [rejection_12_btsp_max_NDM_4_p2_first, T_12_max_p2 > critical_value_12_btsp_max_NDM_4_p2_first];
rejection_12_btsp_sum_NDM_4_p2_first     = [rejection_12_btsp_sum_NDM_4_p2_first, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_4_p2_first];

rejection_12_btsp_max_NDM_5_p2_first     = [rejection_12_btsp_max_NDM_5_p2_first, T_12_max_p2 > critical_value_12_btsp_max_NDM_5_p2_first];
rejection_12_btsp_sum_NDM_5_p2_first     = [rejection_12_btsp_sum_NDM_5_p2_first, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_5_p2_first];

rejection_12_btsp_max_NDM_1_p2_second    = [rejection_12_btsp_max_NDM_1_p2_second, T_12_max_p2 > critical_value_12_btsp_max_NDM_1_p2_second];
rejection_12_btsp_sum_NDM_1_p2_second    = [rejection_12_btsp_sum_NDM_1_p2_second, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_1_p2_second];

rejection_12_btsp_max_NDM_2_p2_second    = [rejection_12_btsp_max_NDM_2_p2_second, T_12_max_p2 > critical_value_12_btsp_max_NDM_2_p2_second];
rejection_12_btsp_sum_NDM_2_p2_second    = [rejection_12_btsp_sum_NDM_2_p2_second, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_2_p2_second];

rejection_12_btsp_max_NDM_3_p2_second    = [rejection_12_btsp_max_NDM_3_p2_second, T_12_max_p2 > critical_value_12_btsp_max_NDM_3_p2_second];
rejection_12_btsp_sum_NDM_3_p2_second    = [rejection_12_btsp_sum_NDM_3_p2_second, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_3_p2_second];

rejection_12_btsp_max_NDM_4_p2_second    = [rejection_12_btsp_max_NDM_4_p2_second, T_12_max_p2 > critical_value_12_btsp_max_NDM_4_p2_second];
rejection_12_btsp_sum_NDM_4_p2_second    = [rejection_12_btsp_sum_NDM_4_p2_second, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_4_p2_second];

rejection_12_btsp_max_NDM_5_p2_second    = [rejection_12_btsp_max_NDM_5_p2_second, T_12_max_p2 > critical_value_12_btsp_max_NDM_5_p2_second];
rejection_12_btsp_sum_NDM_5_p2_second    = [rejection_12_btsp_sum_NDM_5_p2_second, T_12_sum_p2 > critical_value_12_btsp_sum_NDM_5_p2_second];

% (2,1)

rejection_21_btsp_max_LFC_p1     = [rejection_21_btsp_max_LFC_p1, T_21_max_p1 > critical_value_21_btsp_max_LFC_p1];
rejection_21_btsp_sum_LFC_p1     = [rejection_21_btsp_sum_LFC_p1, T_21_sum_p1 > critical_value_21_btsp_sum_LFC_p1];

rejection_21_max_contact_1_p1 =  [rejection_21_max_contact_1_p1, T_21_max_p1 > critical_value_21_max_contact_1_p1];
rejection_21_sum_contact_1_p1 =  [rejection_21_sum_contact_1_p1, T_21_sum_p1 > critical_value_21_sum_contact_1_p1];

rejection_21_max_contact_2_p1 =  [rejection_21_max_contact_2_p1, T_21_max_p1 > critical_value_21_max_contact_2_p1];
rejection_21_sum_contact_2_p1 =  [rejection_21_sum_contact_2_p1, T_21_sum_p1 > critical_value_21_sum_contact_2_p1];

rejection_21_max_contact_3_p1 =  [rejection_21_max_contact_3_p1, T_21_max_p1 > critical_value_21_max_contact_3_p1];
rejection_21_sum_contact_3_p1 =  [rejection_21_sum_contact_3_p1, T_21_sum_p1 > critical_value_21_sum_contact_3_p1];

rejection_21_max_contact_4_p1 =  [rejection_21_max_contact_4_p1, T_21_max_p1 > critical_value_21_max_contact_4_p1];
rejection_21_sum_contact_4_p1 =  [rejection_21_sum_contact_4_p1, T_21_sum_p1 > critical_value_21_sum_contact_4_p1];

rejection_21_max_contact_5_p1 =  [rejection_21_max_contact_5_p1, T_21_max_p1 > critical_value_21_max_contact_5_p1];
rejection_21_sum_contact_5_p1 =  [rejection_21_sum_contact_5_p1, T_21_sum_p1 > critical_value_21_sum_contact_5_p1];

rejection_21_btsp_max_NDM_1_p1     = [rejection_21_btsp_max_NDM_1_p1, T_21_max_p1 > critical_value_21_btsp_max_NDM_1_p1];
rejection_21_btsp_sum_NDM_1_p1     = [rejection_21_btsp_sum_NDM_1_p1, T_21_sum_p1 > critical_value_21_btsp_sum_NDM_1_p1];

rejection_21_btsp_max_NDM_2_p1     = [rejection_21_btsp_max_NDM_2_p1, T_21_max_p1 > critical_value_21_btsp_max_NDM_2_p1];
rejection_21_btsp_sum_NDM_2_p1     = [rejection_21_btsp_sum_NDM_2_p1, T_21_sum_p1 > critical_value_21_btsp_sum_NDM_2_p1];

rejection_21_btsp_max_NDM_3_p1     = [rejection_21_btsp_max_NDM_3_p1, T_21_max_p1 > critical_value_21_btsp_max_NDM_3_p1];
rejection_21_btsp_sum_NDM_3_p1     = [rejection_21_btsp_sum_NDM_3_p1, T_21_sum_p1 > critical_value_21_btsp_sum_NDM_3_p1];

rejection_21_btsp_max_NDM_4_p1     = [rejection_21_btsp_max_NDM_4_p1, T_21_max_p1 > critical_value_21_btsp_max_NDM_4_p1];
rejection_21_btsp_sum_NDM_4_p1     = [rejection_21_btsp_sum_NDM_4_p1, T_21_sum_p1 > critical_value_21_btsp_sum_NDM_4_p1];

rejection_21_btsp_max_NDM_5_p1     = [rejection_21_btsp_max_NDM_5_p1, T_21_max_p1 > critical_value_21_btsp_max_NDM_5_p1];
rejection_21_btsp_sum_NDM_5_p1     = [rejection_21_btsp_sum_NDM_5_p1, T_21_sum_p1 > critical_value_21_btsp_sum_NDM_5_p1];

% p = 2

rejection_21_btsp_max_LFC_p2     = [rejection_21_btsp_max_LFC_p2, T_21_max_p2 > critical_value_21_btsp_max_LFC_p2];
rejection_21_btsp_sum_LFC_p2     = [rejection_21_btsp_sum_LFC_p2, T_21_sum_p2 > critical_value_21_btsp_sum_LFC_p2];

rejection_21_max_contact_1_p2    = [rejection_21_max_contact_1_p2, T_21_max_p2 > critical_value_21_max_contact_1_p2];
rejection_21_sum_contact_1_p2    = [rejection_21_sum_contact_1_p2, T_21_sum_p2 > critical_value_21_sum_contact_1_p2];

rejection_21_max_contact_2_p2    = [rejection_21_max_contact_2_p2, T_21_max_p2 > critical_value_21_max_contact_2_p2];
rejection_21_sum_contact_2_p2    = [rejection_21_sum_contact_2_p2, T_21_sum_p2 > critical_value_21_sum_contact_2_p2];

rejection_21_max_contact_3_p2    = [rejection_21_max_contact_3_p2, T_21_max_p2 > critical_value_21_max_contact_3_p2];
rejection_21_sum_contact_3_p2    = [rejection_21_sum_contact_3_p2, T_21_sum_p2 > critical_value_21_sum_contact_3_p2];

rejection_21_max_contact_4_p2    = [rejection_21_max_contact_4_p2, T_21_max_p2 > critical_value_21_max_contact_4_p2];
rejection_21_sum_contact_4_p2    = [rejection_21_sum_contact_4_p2, T_21_sum_p2 > critical_value_21_sum_contact_4_p2];

rejection_21_max_contact_5_p2    = [rejection_21_max_contact_5_p2, T_21_max_p2 > critical_value_21_max_contact_5_p2];
rejection_21_sum_contact_5_p2    = [rejection_21_sum_contact_5_p2, T_21_sum_p2 > critical_value_21_sum_contact_5_p2];

rejection_21_btsp_max_NDM_1_p2_first     = [rejection_21_btsp_max_NDM_1_p2_first, T_21_max_p2 > critical_value_21_btsp_max_NDM_1_p2_first];
rejection_21_btsp_sum_NDM_1_p2_first     = [rejection_21_btsp_sum_NDM_1_p2_first, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_1_p2_first];

rejection_21_btsp_max_NDM_2_p2_first     = [rejection_21_btsp_max_NDM_2_p2_first, T_21_max_p2 > critical_value_21_btsp_max_NDM_2_p2_first];
rejection_21_btsp_sum_NDM_2_p2_first     = [rejection_21_btsp_sum_NDM_2_p2_first, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_2_p2_first];

rejection_21_btsp_max_NDM_3_p2_first     = [rejection_21_btsp_max_NDM_3_p2_first, T_21_max_p2 > critical_value_21_btsp_max_NDM_3_p2_first];
rejection_21_btsp_sum_NDM_3_p2_first     = [rejection_21_btsp_sum_NDM_3_p2_first, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_3_p2_first];

rejection_21_btsp_max_NDM_4_p2_first     = [rejection_21_btsp_max_NDM_4_p2_first, T_21_max_p2 > critical_value_21_btsp_max_NDM_4_p2_first];
rejection_21_btsp_sum_NDM_4_p2_first     = [rejection_21_btsp_sum_NDM_4_p2_first, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_4_p2_first];

rejection_21_btsp_max_NDM_5_p2_first     = [rejection_21_btsp_max_NDM_5_p2_first, T_21_max_p2 > critical_value_21_btsp_max_NDM_5_p2_first];
rejection_21_btsp_sum_NDM_5_p2_first     = [rejection_21_btsp_sum_NDM_5_p2_first, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_5_p2_first];

rejection_21_btsp_max_NDM_1_p2_second    = [rejection_21_btsp_max_NDM_1_p2_second, T_21_max_p2 > critical_value_21_btsp_max_NDM_1_p2_second];
rejection_21_btsp_sum_NDM_1_p2_second    = [rejection_21_btsp_sum_NDM_1_p2_second, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_1_p2_second];

rejection_21_btsp_max_NDM_2_p2_second    = [rejection_21_btsp_max_NDM_2_p2_second, T_21_max_p2 > critical_value_21_btsp_max_NDM_2_p2_second];
rejection_21_btsp_sum_NDM_2_p2_second    = [rejection_21_btsp_sum_NDM_2_p2_second, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_2_p2_second];

rejection_21_btsp_max_NDM_3_p2_second    = [rejection_21_btsp_max_NDM_3_p2_second, T_21_max_p2 > critical_value_21_btsp_max_NDM_3_p2_second];
rejection_21_btsp_sum_NDM_3_p2_second    = [rejection_21_btsp_sum_NDM_3_p2_second, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_3_p2_second];

rejection_21_btsp_max_NDM_4_p2_second    = [rejection_21_btsp_max_NDM_4_p2_second, T_21_max_p2 > critical_value_21_btsp_max_NDM_4_p2_second];
rejection_21_btsp_sum_NDM_4_p2_second    = [rejection_21_btsp_sum_NDM_4_p2_second, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_4_p2_second];

rejection_21_btsp_max_NDM_5_p2_second    = [rejection_21_btsp_max_NDM_5_p2_second, T_21_max_p2 > critical_value_21_btsp_max_NDM_5_p2_second];
rejection_21_btsp_sum_NDM_5_p2_second    = [rejection_21_btsp_sum_NDM_5_p2_second, T_21_sum_p2 > critical_value_21_btsp_sum_NDM_5_p2_second];

% (2,2)


rejection_22_btsp_max_LFC_p1     = [rejection_22_btsp_max_LFC_p1, T_22_max_p1 > critical_value_22_btsp_max_LFC_p1];
rejection_22_btsp_sum_LFC_p1     = [rejection_22_btsp_sum_LFC_p1, T_22_sum_p1 > critical_value_22_btsp_sum_LFC_p1];

rejection_22_max_contact_1_p1 =  [rejection_22_max_contact_1_p1, T_22_max_p1 > critical_value_22_max_contact_1_p1];
rejection_22_sum_contact_1_p1 =  [rejection_22_sum_contact_1_p1, T_22_sum_p1 > critical_value_22_sum_contact_1_p1];

rejection_22_max_contact_2_p1 =  [rejection_22_max_contact_2_p1, T_22_max_p1 > critical_value_22_max_contact_2_p1];
rejection_22_sum_contact_2_p1 =  [rejection_22_sum_contact_2_p1, T_22_sum_p1 > critical_value_22_sum_contact_2_p1];

rejection_22_max_contact_3_p1 =  [rejection_22_max_contact_3_p1, T_22_max_p1 > critical_value_22_max_contact_3_p1];
rejection_22_sum_contact_3_p1 =  [rejection_22_sum_contact_3_p1, T_22_sum_p1 > critical_value_22_sum_contact_3_p1];

rejection_22_max_contact_4_p1 =  [rejection_22_max_contact_4_p1, T_22_max_p1 > critical_value_22_max_contact_4_p1];
rejection_22_sum_contact_4_p1 =  [rejection_22_sum_contact_4_p1, T_22_sum_p1 > critical_value_22_sum_contact_4_p1];

rejection_22_max_contact_5_p1 =  [rejection_22_max_contact_5_p1, T_22_max_p1 > critical_value_22_max_contact_5_p1];
rejection_22_sum_contact_5_p1 =  [rejection_22_sum_contact_5_p1, T_22_sum_p1 > critical_value_22_sum_contact_5_p1];

rejection_22_btsp_max_NDM_1_p1     = [rejection_22_btsp_max_NDM_1_p1, T_22_max_p1 > critical_value_22_btsp_max_NDM_1_p1];
rejection_22_btsp_sum_NDM_1_p1     = [rejection_22_btsp_sum_NDM_1_p1, T_22_sum_p1 > critical_value_22_btsp_sum_NDM_1_p1];

rejection_22_btsp_max_NDM_2_p1     = [rejection_22_btsp_max_NDM_2_p1, T_22_max_p1 > critical_value_22_btsp_max_NDM_2_p1];
rejection_22_btsp_sum_NDM_2_p1     = [rejection_22_btsp_sum_NDM_2_p1, T_22_sum_p1 > critical_value_22_btsp_sum_NDM_2_p1];

rejection_22_btsp_max_NDM_3_p1     = [rejection_22_btsp_max_NDM_3_p1, T_22_max_p1 > critical_value_22_btsp_max_NDM_3_p1];
rejection_22_btsp_sum_NDM_3_p1     = [rejection_22_btsp_sum_NDM_3_p1, T_22_sum_p1 > critical_value_22_btsp_sum_NDM_3_p1];

rejection_22_btsp_max_NDM_4_p1     = [rejection_22_btsp_max_NDM_4_p1, T_22_max_p1 > critical_value_22_btsp_max_NDM_4_p1];
rejection_22_btsp_sum_NDM_4_p1     = [rejection_22_btsp_sum_NDM_4_p1, T_22_sum_p1 > critical_value_22_btsp_sum_NDM_4_p1];

rejection_22_btsp_max_NDM_5_p1     = [rejection_22_btsp_max_NDM_5_p1, T_22_max_p1 > critical_value_22_btsp_max_NDM_5_p1];
rejection_22_btsp_sum_NDM_5_p1     = [rejection_22_btsp_sum_NDM_5_p1, T_22_sum_p1 > critical_value_22_btsp_sum_NDM_5_p1];

% p = 2

rejection_22_btsp_max_LFC_p2     = [rejection_22_btsp_max_LFC_p2, T_22_max_p2 > critical_value_22_btsp_max_LFC_p2];
rejection_22_btsp_sum_LFC_p2     = [rejection_22_btsp_sum_LFC_p2, T_22_sum_p2 > critical_value_22_btsp_sum_LFC_p2];

rejection_22_max_contact_1_p2    = [rejection_22_max_contact_1_p2, T_22_max_p2 > critical_value_22_max_contact_1_p2];
rejection_22_sum_contact_1_p2    = [rejection_22_sum_contact_1_p2, T_22_sum_p2 > critical_value_22_sum_contact_1_p2];

rejection_22_max_contact_2_p2    = [rejection_22_max_contact_2_p2, T_22_max_p2 > critical_value_22_max_contact_2_p2];
rejection_22_sum_contact_2_p2    = [rejection_22_sum_contact_2_p2, T_22_sum_p2 > critical_value_22_sum_contact_2_p2];

rejection_22_max_contact_3_p2    = [rejection_22_max_contact_3_p2, T_22_max_p2 > critical_value_22_max_contact_3_p2];
rejection_22_sum_contact_3_p2    = [rejection_22_sum_contact_3_p2, T_22_sum_p2 > critical_value_22_sum_contact_3_p2];

rejection_22_max_contact_4_p2    = [rejection_22_max_contact_4_p2, T_22_max_p2 > critical_value_22_max_contact_4_p2];
rejection_22_sum_contact_4_p2    = [rejection_22_sum_contact_4_p2, T_22_sum_p2 > critical_value_22_sum_contact_4_p2];

rejection_22_max_contact_5_p2    = [rejection_22_max_contact_5_p2, T_22_max_p2 > critical_value_22_max_contact_5_p2];
rejection_22_sum_contact_5_p2    = [rejection_22_sum_contact_5_p2, T_22_sum_p2 > critical_value_22_sum_contact_5_p2];

rejection_22_btsp_max_NDM_1_p2_first     = [rejection_22_btsp_max_NDM_1_p2_first, T_22_max_p2 > critical_value_22_btsp_max_NDM_1_p2_first];
rejection_22_btsp_sum_NDM_1_p2_first     = [rejection_22_btsp_sum_NDM_1_p2_first, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_1_p2_first];

rejection_22_btsp_max_NDM_2_p2_first     = [rejection_22_btsp_max_NDM_2_p2_first, T_22_max_p2 > critical_value_22_btsp_max_NDM_2_p2_first];
rejection_22_btsp_sum_NDM_2_p2_first     = [rejection_22_btsp_sum_NDM_2_p2_first, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_2_p2_first];

rejection_22_btsp_max_NDM_3_p2_first     = [rejection_22_btsp_max_NDM_3_p2_first, T_22_max_p2 > critical_value_22_btsp_max_NDM_3_p2_first];
rejection_22_btsp_sum_NDM_3_p2_first     = [rejection_22_btsp_sum_NDM_3_p2_first, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_3_p2_first];

rejection_22_btsp_max_NDM_4_p2_first     = [rejection_22_btsp_max_NDM_4_p2_first, T_22_max_p2 > critical_value_22_btsp_max_NDM_4_p2_first];
rejection_22_btsp_sum_NDM_4_p2_first     = [rejection_22_btsp_sum_NDM_4_p2_first, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_4_p2_first];

rejection_22_btsp_max_NDM_5_p2_first     = [rejection_22_btsp_max_NDM_5_p2_first, T_22_max_p2 > critical_value_22_btsp_max_NDM_5_p2_first];
rejection_22_btsp_sum_NDM_5_p2_first     = [rejection_22_btsp_sum_NDM_5_p2_first, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_5_p2_first];

rejection_22_btsp_max_NDM_1_p2_second    = [rejection_22_btsp_max_NDM_1_p2_second, T_22_max_p2 > critical_value_22_btsp_max_NDM_1_p2_second];
rejection_22_btsp_sum_NDM_1_p2_second    = [rejection_22_btsp_sum_NDM_1_p2_second, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_1_p2_second];

rejection_22_btsp_max_NDM_2_p2_second    = [rejection_22_btsp_max_NDM_2_p2_second, T_22_max_p2 > critical_value_22_btsp_max_NDM_2_p2_second];
rejection_22_btsp_sum_NDM_2_p2_second    = [rejection_22_btsp_sum_NDM_2_p2_second, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_2_p2_second];

rejection_22_btsp_max_NDM_3_p2_second    = [rejection_22_btsp_max_NDM_3_p2_second, T_22_max_p2 > critical_value_22_btsp_max_NDM_3_p2_second];
rejection_22_btsp_sum_NDM_3_p2_second    = [rejection_22_btsp_sum_NDM_3_p2_second, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_3_p2_second];

rejection_22_btsp_max_NDM_4_p2_second    = [rejection_22_btsp_max_NDM_4_p2_second, T_22_max_p2 > critical_value_22_btsp_max_NDM_4_p2_second];
rejection_22_btsp_sum_NDM_4_p2_second    = [rejection_22_btsp_sum_NDM_4_p2_second, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_4_p2_second];

rejection_22_btsp_max_NDM_5_p2_second    = [rejection_22_btsp_max_NDM_5_p2_second, T_22_max_p2 > critical_value_22_btsp_max_NDM_5_p2_second];
rejection_22_btsp_sum_NDM_5_p2_second    = [rejection_22_btsp_sum_NDM_5_p2_second, T_22_sum_p2 > critical_value_22_btsp_sum_NDM_5_p2_second];

% Simulation Ends
end

toc

% % (1,1)

% p = 1

rejection_ratio_11_btsp_max_LFC_p1    = mean(rejection_11_btsp_max_LFC_p1);
rejection_ratio_11_btsp_sum_LFC_p1    = mean(rejection_11_btsp_sum_LFC_p1);

rejection_ratio_11_max_contact_1_p1 = mean(rejection_11_max_contact_1_p1);
rejection_ratio_11_sum_contact_1_p1 = mean(rejection_11_sum_contact_1_p1);

rejection_ratio_11_max_contact_2_p1 = mean(rejection_11_max_contact_2_p1);
rejection_ratio_11_sum_contact_2_p1 = mean(rejection_11_sum_contact_2_p1);

rejection_ratio_11_max_contact_3_p1 = mean(rejection_11_max_contact_3_p1);
rejection_ratio_11_sum_contact_3_p1 = mean(rejection_11_sum_contact_3_p1);

rejection_ratio_11_max_contact_4_p1 = mean(rejection_11_max_contact_4_p1);
rejection_ratio_11_sum_contact_4_p1 = mean(rejection_11_sum_contact_4_p1);

rejection_ratio_11_max_contact_5_p1 = mean(rejection_11_max_contact_5_p1);
rejection_ratio_11_sum_contact_5_p1 = mean(rejection_11_sum_contact_5_p1);

rejection_ratio_11_btsp_max_NDM_1_p1 = mean(rejection_11_btsp_max_NDM_1_p1);
rejection_ratio_11_btsp_sum_NDM_1_p1 = mean(rejection_11_btsp_sum_NDM_1_p1);

rejection_ratio_11_btsp_max_NDM_2_p1 = mean(rejection_11_btsp_max_NDM_2_p1);
rejection_ratio_11_btsp_sum_NDM_2_p1 = mean(rejection_11_btsp_sum_NDM_2_p1);

rejection_ratio_11_btsp_max_NDM_3_p1 = mean(rejection_11_btsp_max_NDM_3_p1);
rejection_ratio_11_btsp_sum_NDM_3_p1 = mean(rejection_11_btsp_sum_NDM_3_p1);

rejection_ratio_11_btsp_max_NDM_4_p1 = mean(rejection_11_btsp_max_NDM_4_p1);
rejection_ratio_11_btsp_sum_NDM_4_p1 = mean(rejection_11_btsp_sum_NDM_4_p1);

rejection_ratio_11_btsp_max_NDM_5_p1 = mean(rejection_11_btsp_max_NDM_5_p1);
rejection_ratio_11_btsp_sum_NDM_5_p1 = mean(rejection_11_btsp_sum_NDM_5_p1);

% p = 2

rejection_ratio_11_btsp_max_LFC_p2   = mean(rejection_11_btsp_max_LFC_p2);
rejection_ratio_11_btsp_sum_LFC_p2   = mean(rejection_11_btsp_sum_LFC_p2);

rejection_ratio_11_max_contact_1_p2 = mean(rejection_11_max_contact_1_p2);
rejection_ratio_11_sum_contact_1_p2 = mean(rejection_11_sum_contact_1_p2);

rejection_ratio_11_max_contact_2_p2 = mean(rejection_11_max_contact_2_p2);
rejection_ratio_11_sum_contact_2_p2 = mean(rejection_11_sum_contact_2_p2);

rejection_ratio_11_max_contact_3_p2 = mean(rejection_11_max_contact_3_p2);
rejection_ratio_11_sum_contact_3_p2 = mean(rejection_11_sum_contact_3_p2);

rejection_ratio_11_max_contact_4_p2 = mean(rejection_11_max_contact_4_p2);
rejection_ratio_11_sum_contact_4_p2 = mean(rejection_11_sum_contact_4_p2);

rejection_ratio_11_max_contact_5_p2 = mean(rejection_11_max_contact_5_p2);
rejection_ratio_11_sum_contact_5_p2 = mean(rejection_11_sum_contact_5_p2);

% first

rejection_ratio_11_btsp_max_NDM_1_p2_first = mean(rejection_11_btsp_max_NDM_1_p2_first);
rejection_ratio_11_btsp_sum_NDM_1_p2_first = mean(rejection_11_btsp_sum_NDM_1_p2_first);

rejection_ratio_11_btsp_max_NDM_2_p2_first = mean(rejection_11_btsp_max_NDM_2_p2_first);
rejection_ratio_11_btsp_sum_NDM_2_p2_first = mean(rejection_11_btsp_sum_NDM_2_p2_first);

rejection_ratio_11_btsp_max_NDM_3_p2_first = mean(rejection_11_btsp_max_NDM_3_p2_first);
rejection_ratio_11_btsp_sum_NDM_3_p2_first = mean(rejection_11_btsp_sum_NDM_3_p2_first);

rejection_ratio_11_btsp_max_NDM_4_p2_first = mean(rejection_11_btsp_max_NDM_4_p2_first);
rejection_ratio_11_btsp_sum_NDM_4_p2_first = mean(rejection_11_btsp_sum_NDM_4_p2_first);

rejection_ratio_11_btsp_max_NDM_5_p2_first = mean(rejection_11_btsp_max_NDM_5_p2_first);
rejection_ratio_11_btsp_sum_NDM_5_p2_first = mean(rejection_11_btsp_sum_NDM_5_p2_first);

% second

rejection_ratio_11_btsp_max_NDM_1_p2_second = mean(rejection_11_btsp_max_NDM_1_p2_second);
rejection_ratio_11_btsp_sum_NDM_1_p2_second = mean(rejection_11_btsp_sum_NDM_1_p2_second);

rejection_ratio_11_btsp_max_NDM_2_p2_second = mean(rejection_11_btsp_max_NDM_2_p2_second);
rejection_ratio_11_btsp_sum_NDM_2_p2_second = mean(rejection_11_btsp_sum_NDM_2_p2_second);

rejection_ratio_11_btsp_max_NDM_3_p2_second = mean(rejection_11_btsp_max_NDM_3_p2_second);
rejection_ratio_11_btsp_sum_NDM_3_p2_second = mean(rejection_11_btsp_sum_NDM_3_p2_second);

rejection_ratio_11_btsp_max_NDM_4_p2_second = mean(rejection_11_btsp_max_NDM_4_p2_second);
rejection_ratio_11_btsp_sum_NDM_4_p2_second = mean(rejection_11_btsp_sum_NDM_4_p2_second);

rejection_ratio_11_btsp_max_NDM_5_p2_second = mean(rejection_11_btsp_max_NDM_5_p2_second);
rejection_ratio_11_btsp_sum_NDM_5_p2_second = mean(rejection_11_btsp_sum_NDM_5_p2_second);

% % (1,2)

% p = 1

rejection_ratio_12_btsp_max_LFC_p1    = mean(rejection_12_btsp_max_LFC_p1);
rejection_ratio_12_btsp_sum_LFC_p1    = mean(rejection_12_btsp_sum_LFC_p1);

rejection_ratio_12_max_contact_1_p1 = mean(rejection_12_max_contact_1_p1);
rejection_ratio_12_sum_contact_1_p1 = mean(rejection_12_sum_contact_1_p1);

rejection_ratio_12_max_contact_2_p1 = mean(rejection_12_max_contact_2_p1);
rejection_ratio_12_sum_contact_2_p1 = mean(rejection_12_sum_contact_2_p1);

rejection_ratio_12_max_contact_3_p1 = mean(rejection_12_max_contact_3_p1);
rejection_ratio_12_sum_contact_3_p1 = mean(rejection_12_sum_contact_3_p1);

rejection_ratio_12_max_contact_4_p1 = mean(rejection_12_max_contact_4_p1);
rejection_ratio_12_sum_contact_4_p1 = mean(rejection_12_sum_contact_4_p1);

rejection_ratio_12_max_contact_5_p1 = mean(rejection_12_max_contact_5_p1);
rejection_ratio_12_sum_contact_5_p1 = mean(rejection_12_sum_contact_5_p1);

rejection_ratio_12_btsp_max_NDM_1_p1 = mean(rejection_12_btsp_max_NDM_1_p1);
rejection_ratio_12_btsp_sum_NDM_1_p1 = mean(rejection_12_btsp_sum_NDM_1_p1);

rejection_ratio_12_btsp_max_NDM_2_p1 = mean(rejection_12_btsp_max_NDM_2_p1);
rejection_ratio_12_btsp_sum_NDM_2_p1 = mean(rejection_12_btsp_sum_NDM_2_p1);

rejection_ratio_12_btsp_max_NDM_3_p1 = mean(rejection_12_btsp_max_NDM_3_p1);
rejection_ratio_12_btsp_sum_NDM_3_p1 = mean(rejection_12_btsp_sum_NDM_3_p1);

rejection_ratio_12_btsp_max_NDM_4_p1 = mean(rejection_12_btsp_max_NDM_4_p1);
rejection_ratio_12_btsp_sum_NDM_4_p1 = mean(rejection_12_btsp_sum_NDM_4_p1);

rejection_ratio_12_btsp_max_NDM_5_p1 = mean(rejection_12_btsp_max_NDM_5_p1);
rejection_ratio_12_btsp_sum_NDM_5_p1 = mean(rejection_12_btsp_sum_NDM_5_p1);


% p = 2

rejection_ratio_12_btsp_max_LFC_p2   = mean(rejection_12_btsp_max_LFC_p2);
rejection_ratio_12_btsp_sum_LFC_p2   = mean(rejection_12_btsp_sum_LFC_p2);

rejection_ratio_12_max_contact_1_p2 = mean(rejection_12_max_contact_1_p2);
rejection_ratio_12_sum_contact_1_p2 = mean(rejection_12_sum_contact_1_p2);

rejection_ratio_12_max_contact_2_p2 = mean(rejection_12_max_contact_2_p2);
rejection_ratio_12_sum_contact_2_p2 = mean(rejection_12_sum_contact_2_p2);

rejection_ratio_12_max_contact_3_p2 = mean(rejection_12_max_contact_3_p2);
rejection_ratio_12_sum_contact_3_p2 = mean(rejection_12_sum_contact_3_p2);

rejection_ratio_12_max_contact_4_p2 = mean(rejection_12_max_contact_4_p2);
rejection_ratio_12_sum_contact_4_p2 = mean(rejection_12_sum_contact_4_p2);

rejection_ratio_12_max_contact_5_p2 = mean(rejection_12_max_contact_5_p2);
rejection_ratio_12_sum_contact_5_p2 = mean(rejection_12_sum_contact_5_p2);

% first

rejection_ratio_12_btsp_max_NDM_1_p2_first = mean(rejection_12_btsp_max_NDM_1_p2_first);
rejection_ratio_12_btsp_sum_NDM_1_p2_first = mean(rejection_12_btsp_sum_NDM_1_p2_first);

rejection_ratio_12_btsp_max_NDM_2_p2_first = mean(rejection_12_btsp_max_NDM_2_p2_first);
rejection_ratio_12_btsp_sum_NDM_2_p2_first = mean(rejection_12_btsp_sum_NDM_2_p2_first);

rejection_ratio_12_btsp_max_NDM_3_p2_first = mean(rejection_12_btsp_max_NDM_3_p2_first);
rejection_ratio_12_btsp_sum_NDM_3_p2_first = mean(rejection_12_btsp_sum_NDM_3_p2_first);

rejection_ratio_12_btsp_max_NDM_4_p2_first = mean(rejection_12_btsp_max_NDM_4_p2_first);
rejection_ratio_12_btsp_sum_NDM_4_p2_first = mean(rejection_12_btsp_sum_NDM_4_p2_first);

rejection_ratio_12_btsp_max_NDM_5_p2_first = mean(rejection_12_btsp_max_NDM_5_p2_first);
rejection_ratio_12_btsp_sum_NDM_5_p2_first = mean(rejection_12_btsp_sum_NDM_5_p2_first);

% second

rejection_ratio_12_btsp_max_NDM_1_p2_second = mean(rejection_12_btsp_max_NDM_1_p2_second);
rejection_ratio_12_btsp_sum_NDM_1_p2_second = mean(rejection_12_btsp_sum_NDM_1_p2_second);

rejection_ratio_12_btsp_max_NDM_2_p2_second = mean(rejection_12_btsp_max_NDM_2_p2_second);
rejection_ratio_12_btsp_sum_NDM_2_p2_second = mean(rejection_12_btsp_sum_NDM_2_p2_second);

rejection_ratio_12_btsp_max_NDM_3_p2_second = mean(rejection_12_btsp_max_NDM_3_p2_second);
rejection_ratio_12_btsp_sum_NDM_3_p2_second = mean(rejection_12_btsp_sum_NDM_3_p2_second);

rejection_ratio_12_btsp_max_NDM_4_p2_second = mean(rejection_12_btsp_max_NDM_4_p2_second);
rejection_ratio_12_btsp_sum_NDM_4_p2_second = mean(rejection_12_btsp_sum_NDM_4_p2_second);

rejection_ratio_12_btsp_max_NDM_5_p2_second = mean(rejection_12_btsp_max_NDM_5_p2_second);
rejection_ratio_12_btsp_sum_NDM_5_p2_second = mean(rejection_12_btsp_sum_NDM_5_p2_second);

% % (2,1)

% p = 1

rejection_ratio_21_btsp_max_LFC_p1    = mean(rejection_21_btsp_max_LFC_p1);
rejection_ratio_21_btsp_sum_LFC_p1    = mean(rejection_21_btsp_sum_LFC_p1);

rejection_ratio_21_max_contact_1_p1 = mean(rejection_21_max_contact_1_p1);
rejection_ratio_21_sum_contact_1_p1 = mean(rejection_21_sum_contact_1_p1);

rejection_ratio_21_max_contact_2_p1 = mean(rejection_21_max_contact_2_p1);
rejection_ratio_21_sum_contact_2_p1 = mean(rejection_21_sum_contact_2_p1);

rejection_ratio_21_max_contact_3_p1 = mean(rejection_21_max_contact_3_p1);
rejection_ratio_21_sum_contact_3_p1 = mean(rejection_21_sum_contact_3_p1);

rejection_ratio_21_max_contact_4_p1 = mean(rejection_21_max_contact_4_p1);
rejection_ratio_21_sum_contact_4_p1 = mean(rejection_21_sum_contact_4_p1);

rejection_ratio_21_max_contact_5_p1 = mean(rejection_21_max_contact_5_p1);
rejection_ratio_21_sum_contact_5_p1 = mean(rejection_21_sum_contact_5_p1);

rejection_ratio_21_btsp_max_NDM_1_p1 = mean(rejection_21_btsp_max_NDM_1_p1);
rejection_ratio_21_btsp_sum_NDM_1_p1 = mean(rejection_21_btsp_sum_NDM_1_p1);

rejection_ratio_21_btsp_max_NDM_2_p1 = mean(rejection_21_btsp_max_NDM_2_p1);
rejection_ratio_21_btsp_sum_NDM_2_p1 = mean(rejection_21_btsp_sum_NDM_2_p1);

rejection_ratio_21_btsp_max_NDM_3_p1 = mean(rejection_21_btsp_max_NDM_3_p1);
rejection_ratio_21_btsp_sum_NDM_3_p1 = mean(rejection_21_btsp_sum_NDM_3_p1);

rejection_ratio_21_btsp_max_NDM_4_p1 = mean(rejection_21_btsp_max_NDM_4_p1);
rejection_ratio_21_btsp_sum_NDM_4_p1 = mean(rejection_21_btsp_sum_NDM_4_p1);

rejection_ratio_21_btsp_max_NDM_5_p1 = mean(rejection_21_btsp_max_NDM_5_p1);
rejection_ratio_21_btsp_sum_NDM_5_p1 = mean(rejection_21_btsp_sum_NDM_5_p1);

% p = 2

rejection_ratio_21_btsp_max_LFC_p2   = mean(rejection_21_btsp_max_LFC_p2);
rejection_ratio_21_btsp_sum_LFC_p2   = mean(rejection_21_btsp_sum_LFC_p2);

rejection_ratio_21_max_contact_1_p2 = mean(rejection_21_max_contact_1_p2);
rejection_ratio_21_sum_contact_1_p2 = mean(rejection_21_sum_contact_1_p2);

rejection_ratio_21_max_contact_2_p2 = mean(rejection_21_max_contact_2_p2);
rejection_ratio_21_sum_contact_2_p2 = mean(rejection_21_sum_contact_2_p2);

rejection_ratio_21_max_contact_3_p2 = mean(rejection_21_max_contact_3_p2);
rejection_ratio_21_sum_contact_3_p2 = mean(rejection_21_sum_contact_3_p2);

rejection_ratio_21_max_contact_4_p2 = mean(rejection_21_max_contact_4_p2);
rejection_ratio_21_sum_contact_4_p2 = mean(rejection_21_sum_contact_4_p2);

rejection_ratio_21_max_contact_5_p2 = mean(rejection_21_max_contact_5_p2);
rejection_ratio_21_sum_contact_5_p2 = mean(rejection_21_sum_contact_5_p2);

% first

rejection_ratio_21_btsp_max_NDM_1_p2_first = mean(rejection_21_btsp_max_NDM_1_p2_first);
rejection_ratio_21_btsp_sum_NDM_1_p2_first = mean(rejection_21_btsp_sum_NDM_1_p2_first);

rejection_ratio_21_btsp_max_NDM_2_p2_first = mean(rejection_21_btsp_max_NDM_2_p2_first);
rejection_ratio_21_btsp_sum_NDM_2_p2_first = mean(rejection_21_btsp_sum_NDM_2_p2_first);

rejection_ratio_21_btsp_max_NDM_3_p2_first = mean(rejection_21_btsp_max_NDM_3_p2_first);
rejection_ratio_21_btsp_sum_NDM_3_p2_first = mean(rejection_21_btsp_sum_NDM_3_p2_first);

rejection_ratio_21_btsp_max_NDM_4_p2_first = mean(rejection_21_btsp_max_NDM_4_p2_first);
rejection_ratio_21_btsp_sum_NDM_4_p2_first = mean(rejection_21_btsp_sum_NDM_4_p2_first);

rejection_ratio_21_btsp_max_NDM_5_p2_first = mean(rejection_21_btsp_max_NDM_5_p2_first);
rejection_ratio_21_btsp_sum_NDM_5_p2_first = mean(rejection_21_btsp_sum_NDM_5_p2_first);

% second

rejection_ratio_21_btsp_max_NDM_1_p2_second = mean(rejection_21_btsp_max_NDM_1_p2_second);
rejection_ratio_21_btsp_sum_NDM_1_p2_second = mean(rejection_21_btsp_sum_NDM_1_p2_second);
% 
rejection_ratio_21_btsp_max_NDM_2_p2_second = mean(rejection_21_btsp_max_NDM_2_p2_second);
rejection_ratio_21_btsp_sum_NDM_2_p2_second = mean(rejection_21_btsp_sum_NDM_2_p2_second);

rejection_ratio_21_btsp_max_NDM_3_p2_second = mean(rejection_21_btsp_max_NDM_3_p2_second);
rejection_ratio_21_btsp_sum_NDM_3_p2_second = mean(rejection_21_btsp_sum_NDM_3_p2_second);

rejection_ratio_21_btsp_max_NDM_4_p2_second = mean(rejection_21_btsp_max_NDM_4_p2_second);
rejection_ratio_21_btsp_sum_NDM_4_p2_second = mean(rejection_21_btsp_sum_NDM_4_p2_second);

rejection_ratio_21_btsp_max_NDM_5_p2_second = mean(rejection_21_btsp_max_NDM_5_p2_second);
rejection_ratio_21_btsp_sum_NDM_5_p2_second = mean(rejection_21_btsp_sum_NDM_5_p2_second);

% % (2,2)

% p = 1

% % (1,1)

% p = 1

rejection_ratio_22_btsp_max_LFC_p1    = mean(rejection_22_btsp_max_LFC_p1);
rejection_ratio_22_btsp_sum_LFC_p1    = mean(rejection_22_btsp_sum_LFC_p1);

rejection_ratio_22_max_contact_1_p1 = mean(rejection_22_max_contact_1_p1);
rejection_ratio_22_sum_contact_1_p1 = mean(rejection_22_sum_contact_1_p1);

rejection_ratio_22_max_contact_2_p1 = mean(rejection_22_max_contact_2_p1);
rejection_ratio_22_sum_contact_2_p1 = mean(rejection_22_sum_contact_2_p1);

rejection_ratio_22_max_contact_3_p1 = mean(rejection_22_max_contact_3_p1);
rejection_ratio_22_sum_contact_3_p1 = mean(rejection_22_sum_contact_3_p1);

rejection_ratio_22_max_contact_4_p1 = mean(rejection_22_max_contact_4_p1);
rejection_ratio_22_sum_contact_4_p1 = mean(rejection_22_sum_contact_4_p1);

rejection_ratio_22_max_contact_5_p1 = mean(rejection_22_max_contact_5_p1);
rejection_ratio_22_sum_contact_5_p1 = mean(rejection_22_sum_contact_5_p1);

rejection_ratio_22_btsp_max_NDM_1_p1 = mean(rejection_22_btsp_max_NDM_1_p1);
rejection_ratio_22_btsp_sum_NDM_1_p1 = mean(rejection_22_btsp_sum_NDM_1_p1);

rejection_ratio_22_btsp_max_NDM_2_p1 = mean(rejection_22_btsp_max_NDM_2_p1);
rejection_ratio_22_btsp_sum_NDM_2_p1 = mean(rejection_22_btsp_sum_NDM_2_p1);

rejection_ratio_22_btsp_max_NDM_3_p1 = mean(rejection_22_btsp_max_NDM_3_p1);
rejection_ratio_22_btsp_sum_NDM_3_p1 = mean(rejection_22_btsp_sum_NDM_3_p1);

rejection_ratio_22_btsp_max_NDM_4_p1 = mean(rejection_22_btsp_max_NDM_4_p1);
rejection_ratio_22_btsp_sum_NDM_4_p1 = mean(rejection_22_btsp_sum_NDM_4_p1);

rejection_ratio_22_btsp_max_NDM_5_p1 = mean(rejection_22_btsp_max_NDM_5_p1);
rejection_ratio_22_btsp_sum_NDM_5_p1 = mean(rejection_22_btsp_sum_NDM_5_p1);

% p = 2

rejection_ratio_22_btsp_max_LFC_p2   = mean(rejection_22_btsp_max_LFC_p2);
rejection_ratio_22_btsp_sum_LFC_p2   = mean(rejection_22_btsp_sum_LFC_p2);

rejection_ratio_22_max_contact_1_p2 = mean(rejection_22_max_contact_1_p2);
rejection_ratio_22_sum_contact_1_p2 = mean(rejection_22_sum_contact_1_p2);

rejection_ratio_22_max_contact_2_p2 = mean(rejection_22_max_contact_2_p2);
rejection_ratio_22_sum_contact_2_p2 = mean(rejection_22_sum_contact_2_p2);

rejection_ratio_22_max_contact_3_p2 = mean(rejection_22_max_contact_3_p2);
rejection_ratio_22_sum_contact_3_p2 = mean(rejection_22_sum_contact_3_p2);

rejection_ratio_22_max_contact_4_p2 = mean(rejection_22_max_contact_4_p2);
rejection_ratio_22_sum_contact_4_p2 = mean(rejection_22_sum_contact_4_p2);

rejection_ratio_22_max_contact_5_p2 = mean(rejection_22_max_contact_5_p2);
rejection_ratio_22_sum_contact_5_p2 = mean(rejection_22_sum_contact_5_p2);

% first

rejection_ratio_22_btsp_max_NDM_1_p2_first = mean(rejection_22_btsp_max_NDM_1_p2_first);
rejection_ratio_22_btsp_sum_NDM_1_p2_first = mean(rejection_22_btsp_sum_NDM_1_p2_first);

rejection_ratio_22_btsp_max_NDM_2_p2_first = mean(rejection_22_btsp_max_NDM_2_p2_first);
rejection_ratio_22_btsp_sum_NDM_2_p2_first = mean(rejection_22_btsp_sum_NDM_2_p2_first);

rejection_ratio_22_btsp_max_NDM_3_p2_first = mean(rejection_22_btsp_max_NDM_3_p2_first);
rejection_ratio_22_btsp_sum_NDM_3_p2_first = mean(rejection_22_btsp_sum_NDM_3_p2_first);

rejection_ratio_22_btsp_max_NDM_4_p2_first = mean(rejection_22_btsp_max_NDM_4_p2_first);
rejection_ratio_22_btsp_sum_NDM_4_p2_first = mean(rejection_22_btsp_sum_NDM_4_p2_first);

rejection_ratio_22_btsp_max_NDM_5_p2_first = mean(rejection_22_btsp_max_NDM_5_p2_first);
rejection_ratio_22_btsp_sum_NDM_5_p2_first = mean(rejection_22_btsp_sum_NDM_5_p2_first);

% second

rejection_ratio_22_btsp_max_NDM_1_p2_second = mean(rejection_22_btsp_max_NDM_1_p2_second);
rejection_ratio_22_btsp_sum_NDM_1_p2_second = mean(rejection_22_btsp_sum_NDM_1_p2_second);

rejection_ratio_22_btsp_max_NDM_2_p2_second = mean(rejection_22_btsp_max_NDM_2_p2_second);
rejection_ratio_22_btsp_sum_NDM_2_p2_second = mean(rejection_22_btsp_sum_NDM_2_p2_second);

rejection_ratio_22_btsp_max_NDM_3_p2_second = mean(rejection_22_btsp_max_NDM_3_p2_second);
rejection_ratio_22_btsp_sum_NDM_3_p2_second = mean(rejection_22_btsp_sum_NDM_3_p2_second);

rejection_ratio_22_btsp_max_NDM_4_p2_second = mean(rejection_22_btsp_max_NDM_4_p2_second);
rejection_ratio_22_btsp_sum_NDM_4_p2_second = mean(rejection_22_btsp_sum_NDM_4_p2_second);

rejection_ratio_22_btsp_max_NDM_5_p2_second = mean(rejection_22_btsp_max_NDM_5_p2_second);
rejection_ratio_22_btsp_sum_NDM_5_p2_second = mean(rejection_22_btsp_sum_NDM_5_p2_second);

Results_max_p1 = [Results_max_p1; 
     
     [ N design "(1,1)" rejection_ratio_11_btsp_max_LFC_p1 ...
     rejection_ratio_11_max_contact_1_p1    rejection_ratio_11_max_contact_2_p1    rejection_ratio_11_max_contact_3_p1    rejection_ratio_11_max_contact_4_p1 rejection_ratio_11_max_contact_5_p1     ... 
     rejection_ratio_11_btsp_max_NDM_1_p1   rejection_ratio_11_btsp_max_NDM_2_p1   rejection_ratio_11_btsp_max_NDM_3_p1   rejection_ratio_11_btsp_max_NDM_4_p1   rejection_ratio_11_btsp_max_NDM_5_p1 ...
     ];     

     [ N design "(1,2)" rejection_ratio_12_btsp_max_LFC_p1 ...
     rejection_ratio_12_max_contact_1_p1    rejection_ratio_12_max_contact_2_p1    rejection_ratio_12_max_contact_3_p1    rejection_ratio_12_max_contact_4_p1 rejection_ratio_12_max_contact_5_p1    ... 
     rejection_ratio_12_btsp_max_NDM_1_p1   rejection_ratio_12_btsp_max_NDM_2_p1   rejection_ratio_12_btsp_max_NDM_3_p1   rejection_ratio_12_btsp_max_NDM_4_p1   rejection_ratio_12_btsp_max_NDM_5_p1 ...
     ];     

     [ N design "(2,1)" rejection_ratio_21_btsp_max_LFC_p1 ...
     rejection_ratio_21_max_contact_1_p1    rejection_ratio_21_max_contact_2_p1    rejection_ratio_21_max_contact_3_p1    rejection_ratio_21_max_contact_4_p1 rejection_ratio_21_max_contact_5_p1   ... 
     rejection_ratio_21_btsp_max_NDM_1_p1   rejection_ratio_21_btsp_max_NDM_2_p1   rejection_ratio_21_btsp_max_NDM_3_p1   rejection_ratio_21_btsp_max_NDM_4_p1   rejection_ratio_21_btsp_max_NDM_5_p1 ...
     ];     

     [ N design "(2,2)" rejection_ratio_22_btsp_max_LFC_p1 ...
     rejection_ratio_22_max_contact_1_p1    rejection_ratio_22_max_contact_2_p1    rejection_ratio_22_max_contact_3_p1    rejection_ratio_22_max_contact_4_p1 rejection_ratio_22_max_contact_5_p1  ... 
     rejection_ratio_22_btsp_max_NDM_1_p1   rejection_ratio_22_btsp_max_NDM_2_p1   rejection_ratio_22_btsp_max_NDM_3_p1   rejection_ratio_22_btsp_max_NDM_4_p1   rejection_ratio_22_btsp_max_NDM_5_p1 ...
     ]  
     ]
 
Results_max_p2 = [ Results_max_p2; 
     
     [ N design "(1,1)" rejection_ratio_11_btsp_max_LFC_p2 ...
     rejection_ratio_11_max_contact_1_p2    rejection_ratio_11_max_contact_2_p2    rejection_ratio_11_max_contact_3_p2    rejection_ratio_11_max_contact_4_p2  rejection_ratio_11_max_contact_5_p2  ... 
     rejection_ratio_11_btsp_max_NDM_1_p2_first  rejection_ratio_11_btsp_max_NDM_2_p2_first  rejection_ratio_11_btsp_max_NDM_3_p2_first  rejection_ratio_11_btsp_max_NDM_4_p2_first  rejection_ratio_11_btsp_max_NDM_5_p2_first ...
     rejection_ratio_11_btsp_max_NDM_1_p2_second rejection_ratio_11_btsp_max_NDM_2_p2_second rejection_ratio_11_btsp_max_NDM_3_p2_second rejection_ratio_11_btsp_max_NDM_4_p2_second rejection_ratio_11_btsp_max_NDM_5_p2_second ...
    ];
       
     [ N design "(1,2)" rejection_ratio_12_btsp_max_LFC_p2 ...
     rejection_ratio_12_max_contact_1_p2    rejection_ratio_12_max_contact_2_p2    rejection_ratio_12_max_contact_3_p2    rejection_ratio_12_max_contact_4_p2 rejection_ratio_12_max_contact_5_p2    ... 
     rejection_ratio_12_btsp_max_NDM_1_p2_first  rejection_ratio_12_btsp_max_NDM_2_p2_first  rejection_ratio_12_btsp_max_NDM_3_p2_first  rejection_ratio_12_btsp_max_NDM_4_p2_first  rejection_ratio_12_btsp_max_NDM_5_p2_first ...
     rejection_ratio_12_btsp_max_NDM_1_p2_second rejection_ratio_12_btsp_max_NDM_2_p2_second rejection_ratio_12_btsp_max_NDM_3_p2_second rejection_ratio_12_btsp_max_NDM_4_p2_second rejection_ratio_12_btsp_max_NDM_5_p2_second ...
     ];

     [ N design "(2,1)" rejection_ratio_21_btsp_max_LFC_p2 ...
     rejection_ratio_21_max_contact_1_p2    rejection_ratio_21_max_contact_2_p2    rejection_ratio_21_max_contact_3_p2    rejection_ratio_21_max_contact_4_p2 rejection_ratio_21_max_contact_5_p2    ... 
     rejection_ratio_21_btsp_max_NDM_1_p2_first  rejection_ratio_21_btsp_max_NDM_2_p2_first  rejection_ratio_21_btsp_max_NDM_3_p2_first  rejection_ratio_21_btsp_max_NDM_4_p2_first  rejection_ratio_21_btsp_max_NDM_5_p2_first ...
     rejection_ratio_21_btsp_max_NDM_1_p2_second rejection_ratio_21_btsp_max_NDM_2_p2_second rejection_ratio_21_btsp_max_NDM_3_p2_second rejection_ratio_21_btsp_max_NDM_4_p2_second rejection_ratio_21_btsp_max_NDM_5_p2_second ...
     ];
 
     [ N design "(2,2)" rejection_ratio_22_btsp_max_LFC_p2 ...
     rejection_ratio_22_max_contact_1_p2    rejection_ratio_22_max_contact_2_p2    rejection_ratio_22_max_contact_3_p2    rejection_ratio_22_max_contact_4_p2 rejection_ratio_22_max_contact_5_p2   ... 
     rejection_ratio_22_btsp_max_NDM_1_p2_first  rejection_ratio_22_btsp_max_NDM_2_p2_first  rejection_ratio_22_btsp_max_NDM_3_p2_first  rejection_ratio_22_btsp_max_NDM_4_p2_first  rejection_ratio_22_btsp_max_NDM_5_p2_first ...
     rejection_ratio_22_btsp_max_NDM_1_p2_second rejection_ratio_22_btsp_max_NDM_2_p2_second rejection_ratio_22_btsp_max_NDM_3_p2_second rejection_ratio_22_btsp_max_NDM_4_p2_second rejection_ratio_22_btsp_max_NDM_5_p2_second ...
     ];      
     ]
 


Results_sum_p1 = [Results_sum_p1; 
     
     [ N design "(1,1)" rejection_ratio_11_btsp_sum_LFC_p1 ...
     rejection_ratio_11_sum_contact_1_p1    rejection_ratio_11_sum_contact_2_p1    rejection_ratio_11_sum_contact_3_p1    rejection_ratio_11_sum_contact_4_p1 rejection_ratio_11_sum_contact_5_p1     ... 
     rejection_ratio_11_btsp_sum_NDM_1_p1   rejection_ratio_11_btsp_sum_NDM_2_p1   rejection_ratio_11_btsp_sum_NDM_3_p1   rejection_ratio_11_btsp_sum_NDM_4_p1   rejection_ratio_11_btsp_sum_NDM_5_p1 ...
     ];     

     [ N design "(1,2)" rejection_ratio_12_btsp_sum_LFC_p1 ...
     rejection_ratio_12_sum_contact_1_p1    rejection_ratio_12_sum_contact_2_p1    rejection_ratio_12_sum_contact_3_p1    rejection_ratio_12_sum_contact_4_p1 rejection_ratio_12_sum_contact_5_p1    ... 
     rejection_ratio_12_btsp_sum_NDM_1_p1   rejection_ratio_12_btsp_sum_NDM_2_p1   rejection_ratio_12_btsp_sum_NDM_3_p1   rejection_ratio_12_btsp_sum_NDM_4_p1   rejection_ratio_12_btsp_sum_NDM_5_p1 ...
     ];     

     [ N design "(2,1)" rejection_ratio_21_btsp_sum_LFC_p1 ...
     rejection_ratio_21_sum_contact_1_p1    rejection_ratio_21_sum_contact_2_p1    rejection_ratio_21_sum_contact_3_p1    rejection_ratio_21_sum_contact_4_p1 rejection_ratio_21_sum_contact_5_p1   ... 
     rejection_ratio_21_btsp_sum_NDM_1_p1   rejection_ratio_21_btsp_sum_NDM_2_p1   rejection_ratio_21_btsp_sum_NDM_3_p1   rejection_ratio_21_btsp_sum_NDM_4_p1   rejection_ratio_21_btsp_sum_NDM_5_p1 ...
     ];     

     [ N design "(2,2)" rejection_ratio_22_btsp_sum_LFC_p1 ...
     rejection_ratio_22_sum_contact_1_p1    rejection_ratio_22_sum_contact_2_p1    rejection_ratio_22_sum_contact_3_p1    rejection_ratio_22_sum_contact_4_p1 rejection_ratio_22_sum_contact_5_p1  ... 
     rejection_ratio_22_btsp_sum_NDM_1_p1   rejection_ratio_22_btsp_sum_NDM_2_p1   rejection_ratio_22_btsp_sum_NDM_3_p1   rejection_ratio_22_btsp_sum_NDM_4_p1   rejection_ratio_22_btsp_sum_NDM_5_p1 ...
     ]  
     ]
 
Results_sum_p2 = [ Results_sum_p2; 
     
     [ N design "(1,1)" rejection_ratio_11_btsp_sum_LFC_p2 ...
     rejection_ratio_11_sum_contact_1_p2    rejection_ratio_11_sum_contact_2_p2    rejection_ratio_11_sum_contact_3_p2    rejection_ratio_11_sum_contact_4_p2  rejection_ratio_11_sum_contact_5_p2  ... 
     rejection_ratio_11_btsp_sum_NDM_1_p2_first  rejection_ratio_11_btsp_sum_NDM_2_p2_first  rejection_ratio_11_btsp_sum_NDM_3_p2_first  rejection_ratio_11_btsp_sum_NDM_4_p2_first  rejection_ratio_11_btsp_sum_NDM_5_p2_first ...
     rejection_ratio_11_btsp_sum_NDM_1_p2_second rejection_ratio_11_btsp_sum_NDM_2_p2_second rejection_ratio_11_btsp_sum_NDM_3_p2_second rejection_ratio_11_btsp_sum_NDM_4_p2_second rejection_ratio_11_btsp_sum_NDM_5_p2_second ...
    ];
       
     [ N design "(1,2)" rejection_ratio_12_btsp_sum_LFC_p2 ...
     rejection_ratio_12_sum_contact_1_p2    rejection_ratio_12_sum_contact_2_p2    rejection_ratio_12_sum_contact_3_p2    rejection_ratio_12_sum_contact_4_p2 rejection_ratio_12_sum_contact_5_p2    ... 
     rejection_ratio_12_btsp_sum_NDM_1_p2_first  rejection_ratio_12_btsp_sum_NDM_2_p2_first  rejection_ratio_12_btsp_sum_NDM_3_p2_first  rejection_ratio_12_btsp_sum_NDM_4_p2_first  rejection_ratio_12_btsp_sum_NDM_5_p2_first ...
     rejection_ratio_12_btsp_sum_NDM_1_p2_second rejection_ratio_12_btsp_sum_NDM_2_p2_second rejection_ratio_12_btsp_sum_NDM_3_p2_second rejection_ratio_12_btsp_sum_NDM_4_p2_second rejection_ratio_12_btsp_sum_NDM_5_p2_second ...
     ];

     [ N design "(2,1)" rejection_ratio_21_btsp_sum_LFC_p2 ...
     rejection_ratio_21_sum_contact_1_p2    rejection_ratio_21_sum_contact_2_p2    rejection_ratio_21_sum_contact_3_p2    rejection_ratio_21_sum_contact_4_p2 rejection_ratio_21_sum_contact_5_p2    ... 
     rejection_ratio_21_btsp_sum_NDM_1_p2_first  rejection_ratio_21_btsp_sum_NDM_2_p2_first  rejection_ratio_21_btsp_sum_NDM_3_p2_first  rejection_ratio_21_btsp_sum_NDM_4_p2_first  rejection_ratio_21_btsp_sum_NDM_5_p2_first ...
     rejection_ratio_21_btsp_sum_NDM_1_p2_second rejection_ratio_21_btsp_sum_NDM_2_p2_second rejection_ratio_21_btsp_sum_NDM_3_p2_second rejection_ratio_21_btsp_sum_NDM_4_p2_second rejection_ratio_21_btsp_sum_NDM_5_p2_second ...
     ];
 
     [ N design "(2,2)" rejection_ratio_22_btsp_sum_LFC_p2 ...
     rejection_ratio_22_sum_contact_1_p2    rejection_ratio_22_sum_contact_2_p2    rejection_ratio_22_sum_contact_3_p2    rejection_ratio_22_sum_contact_4_p2 rejection_ratio_22_sum_contact_5_p2   ... 
     rejection_ratio_22_btsp_sum_NDM_1_p2_first  rejection_ratio_22_btsp_sum_NDM_2_p2_first  rejection_ratio_22_btsp_sum_NDM_3_p2_first  rejection_ratio_22_btsp_sum_NDM_4_p2_first  rejection_ratio_22_btsp_sum_NDM_5_p2_first ...
     rejection_ratio_22_btsp_sum_NDM_1_p2_second rejection_ratio_22_btsp_sum_NDM_2_p2_second rejection_ratio_22_btsp_sum_NDM_3_p2_second rejection_ratio_22_btsp_sum_NDM_4_p2_second rejection_ratio_22_btsp_sum_NDM_5_p2_second ...
     ];      
     ]
 

writematrix(Results_max_p1,'Simul_Unif_p1_max.csv')
writematrix(Results_sum_p1,'Simul_Unif_p1_sum.csv')
writematrix(Results_max_p2,'Simul_Unif_p2_max.csv')
writematrix(Results_sum_p2,'Simul_Unif_p2_sum.csv')

toc

end
end
