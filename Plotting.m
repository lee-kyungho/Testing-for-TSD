
% Simulation Arguments

N = 100000;
n_simul = 1;
T = 4;
ngrid = 100;
p = 1;

path = "Figures/crossing_normal/";

% Sampling

for design = 1:12

X = zeros(N,T+1,n_simul);
X0 = randn(N,1,n_simul);
X(:,1,:) = X0;

for i = 1:T
    X(:,i+1,:) = 0.5*X(:,i,:) +  randn(N,1,n_simul);
end

Y = zeros(N,T+1,n_simul);
Y0 = randn(N,1,n_simul);
Y(:,1,:) = Y0;

for i = 1:T
    Y(:,i+1,:) = 0.5*Y(:,i,:) + randn(N,1,n_simul);
end

if design == 1
Y(:,1,:) = Y(:,1,:)*0.5;
elseif design == 2
Y(:,1,:) = Y(:,1,:)*2;
elseif design == 3
Y(:,3,:) = Y(:,3,:)*0.5;
elseif design == 4
Y(:,3,:) = Y(:,3,:)*2;
elseif design == 5
Y(:,T+1,:) = Y(:,T+1,:)*0.5;
elseif design == 6
Y(:,T+1,:) = Y(:,T+1,:)*2;
elseif design == 7
Y(:,T+1,:) = Y(:,T+1,:)*0.25;
elseif design == 8
Y(:,T+1,:) = Y(:,T+1,:)*4;
elseif design == 9
Y(:,T+1,:) = Y(:,T+1,:)*0.25;
elseif design == 10
Y(:,T+1,:) = Y(:,T+1,:)*4;
elseif design == 11
Y(:,T+1,:) = Y(:,T+1,:)*0.25;
elseif design == 12
Y(:,T+1,:) = Y(:,T+1,:)*4;

end

r_N = sqrt(size(X,1)*size(Y,1)/(size(X,1) + size(Y,1)));
% r_N = sqrt(size(X,1));

R = 1;
    
sample1 = squeeze(X(:,:,R));
sample2 = squeeze(Y(:,:,R));

grid = linspace(min(min(sample1,sample2),[],'all'),max(max(sample1,sample2),[],'all'),ngrid)';

% Calculate Test Statistics

% Note on dimension on D

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

Lamb_11_max = Lambda(D_11,p,'max'); % Output dim: ngrid * 1
Lamb_11_sum = Lambda(D_11,p,'sum'); % Output dim: ngrid * 1

Lamb_12_max = Lambda(D_12,p,'max'); % Output dim: ngrid * 1
Lamb_12_sum = Lambda(D_12,p,'sum'); % Output dim: ngrid * 1

Lamb_21_max = Lambda(D_21,p,'max'); % Output dim: ngrid * 1
Lamb_21_sum = Lambda(D_21,p,'sum'); % Output dim: ngrid * 1

% For n=2, m=2, we need to make collection of D_22 and D_12_T

D_22_collection = cat(2,D_22,D_12_T); % ngrid * J

Lamb_22_max = Lambda(D_22_collection,p,'max'); % Output dim: ngrid * 1
Lamb_22_sum = Lambda(D_22_collection,p,'sum'); % Output dim: ngrid * 1

%---- Integration for T_{N} ----------%

T_11_max = r_N^p * trapz(Lamb_11_max);
T_11_sum = r_N^p * trapz(Lamb_11_sum);

T_12_max = r_N^p * trapz(Lamb_12_max);
T_12_sum = r_N^p * trapz(Lamb_12_sum);

T_21_max = r_N^p * trapz(Lamb_21_max);
T_21_sum = r_N^p * trapz(Lamb_21_sum);

T_22_max = r_N^p * trapz(Lamb_22_max);
T_22_sum = r_N^p * trapz(Lamb_22_sum);


for i = 1:size(D_11,2)
plot(grid,D_11(:,i));
hold on
end
xlabel("x")
ylabel("D^{(1,1)}")
legend("t=0","t=1","t=2","t=3","t=4")
saveas(gcf, path + "Design, " + design + ' D_11_DGP_'+ date +'.png')
clf('reset')

for i = 1:size(D_12,2)
plot(grid,D_12(:,i));
hold on
end
xlabel("x")
ylabel("D^{(1,2)}")
legend("t=0","t=1","t=2","t=3","t=4")
saveas(gcf, path + "Design, " + design + ' D_12_DGP_'+ date +'.png')
clf('reset')


for i = 1:size(D_21,2)
plot(grid,D_21(:,i));
hold on
end
xlabel("x")
ylabel("D^{(2,1)}")
legend("t=0","t=1","t=2","t=3","t=4")
saveas(gcf, path + "Design, " + design + ' D_21_DGP_'+ date +'.png')
clf('reset')

for i = 1:size(D_22_collection,2)
plot(grid,D_22_collection(:,i));
hold on
end
xlabel("x")
ylabel("D^{(2,2)}(x,t) / D^{(1,2)}(x,T) ")
legend("t=0","t=1","t=2","t=3","t=4","Terminal")
saveas(gcf, path + "Design, " + design + ' D_22_collection_DGP_'+ date +'.png')
clf('reset')


plot(grid,Lamb_11_sum)
xlabel("x")
ylabel("\Lambda^{(1,1)} Sum")
saveas(gcf, path + "Design, " + design + ' Lamb_11_sum_DGP_'+ date +'.png')
clf('reset')

plot(grid,Lamb_12_sum)
xlabel("x")
ylabel("\Lambda^{(1,2)} Sum")
saveas(gcf, path + "Design, " + design + ' Lamb_12_sum_DGP_'+ date +'.png')
clf('reset')

plot(grid,Lamb_21_sum)
xlabel("x")
ylabel("\Lambda^{(2,1)} Sum")
saveas(gcf, path + "Design, " + design + ' Lamb_21_sum_DGP_'+ date +'.png')
clf('reset')

plot(grid,Lamb_22_sum)
xlabel("x")
ylabel("\Lambda^{(2,2)} Sum")
saveas(gcf, path + "Design, " + design + ' Lamb_22_sum_DGP_'+ date +'.png')
clf('reset')


plot(grid,Lamb_11_sum)
xlabel("x")
ylabel("\Lambda^{(1,1)} Sum")
saveas(gcf, path + "Design, " + design + ' Lamb_11_sum_DGP_'+ date +'.png')
clf('reset')

plot(grid,Lamb_12_sum)
xlabel("x")
ylabel("\Lambda^{(1,2)} Sum")
saveas(gcf, path + "Design, " + design + ' Lamb_12_sum_DGP_'+ date +'.png')
clf('reset')

plot(grid,Lamb_21_sum)
xlabel("x")
ylabel("\Lambda^{(2,1)} Sum")
saveas(gcf, path + "Design, " + design + ' Lamb_21_sum_DGP_'+ date +'.png')
clf('reset')

plot(grid,Lamb_22_sum)
xlabel("x")
ylabel("\Lambda^{(2,2)} Sum")
saveas(gcf, path + "Design, " + design + ' Lamb_22_sum_DGP_'+ date +'.png')
clf('reset')


i = 1;
for n = 1:2
    for m = 1:2
    i = i+1;
operated_X = squeeze(operation(n,m,X,grid));
operated_Y = squeeze(operation(n,m,Y,grid));

set(plot(grid,operated_X(:,1)), 'LineStyle', '-.', 'Color', [0,0.0,0.5],'LineWidth',1);
hold on
set(plot(grid,operated_X(:,2)), 'LineStyle', '-.', 'Color', [0,0.1,0.9],'LineWidth',1);
hold on
set(plot(grid,operated_X(:,3)), 'LineStyle', '-.', 'Color', [0,0.3,0.7],'LineWidth',1);
hold on
set(plot(grid,operated_X(:,4)), 'LineStyle', '-.', 'Color', [0,0.3,0.7],'LineWidth',1);
hold on
set(plot(grid,operated_X(:,5)), 'LineStyle', '-.', 'Color', [0,0.3,0.7],'LineWidth',1);
% hold on
set(plot(grid,operated_Y(:,1)), 'LineStyle', '-', 'Color', [0.5,0,0],'LineWidth',1);
hold on
set(plot(grid,operated_Y(:,2)), 'LineStyle', '-', 'Color', [0.8,0.2,0],'LineWidth',1);
hold on
set(plot(grid,operated_Y(:,3)), 'LineStyle', '-', 'Color', [0.7,0.1,0.2],'LineWidth',1);
hold on
set(plot(grid,operated_Y(:,4)), 'LineStyle', '-', 'Color', [0.7,0.2,0],'LineWidth',1);
hold on
set(plot(grid,operated_Y(:,5)), 'LineStyle', '-', 'Color', [0.8,0.2,0],'LineWidth',1);
xlabel("z")
ylabel("Operator(z)")
legend('X1,0','X1,1','X1,2','X1,3','X1,4','X2,0','X2,1','X2,2','X2,3','X2,4','Location','northwest');
title("Order = (" + n + "," + m + ")");
saveas(gcf, path + "Design, " + design +  " Order = (" + n + "," + m + ")" + date +'.png')
clf('reset')
    end
end
end
