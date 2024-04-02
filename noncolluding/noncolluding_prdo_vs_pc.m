%% D2D中断概率 vs P_c
% case 1
clear;
lambda = 1;      % 信道系数方差
sigma2 = -10;    % dBm，噪声功率
P_c_array = 0:2:40;         % dBm
P_dmax = 30;     % dBm
R_d = 0.1;       % Mbps     % case 1

% 理论值
PR_do_theo_array = zeros(1,numel(P_c_array));
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
    A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    PR_do_theo_array(idx) = PR_do;
end
plot(P_c_array, PR_do_theo_array, 'r.-', 'LineWidth',1.0);
hold on;

% 仿真值
PR_do_simu_array = zeros(1,numel(P_c_array));
M = 10000;
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    sum_do = 0;     % D2D中断总次数
    for m = 1:M
        % 生成服从均匀分布的随机功率
        P_d = dbm2w(P_dmax) * rand();
        % 生成信道
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % 计算SINR
        SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
        % 计算速率
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
    end
    PR_do_simu_array(idx) = sum_do / M;
end
plot(P_c_array, PR_do_simu_array, 'rs', 'LineWidth',1.0);
hold on;


% case 2
R_d = 0.5;       % Mbps

% 理论值
PR_do_theo_array = zeros(1,numel(P_c_array));
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
    A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    PR_do_theo_array(idx) = PR_do;
end
plot(P_c_array, PR_do_theo_array, 'g.-', 'LineWidth',1.0);
hold on;

% 仿真值
PR_do_simu_array = zeros(1,numel(P_c_array));
M = 10000;
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    sum_do = 0;     % D2D中断总次数
    for m = 1:M
        % 生成服从均匀分布的随机功率
        P_d = dbm2w(P_dmax) * rand();
        % 生成信道
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % 计算SINR
        SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
        % 计算速率
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
    end
    PR_do_simu_array(idx) = sum_do / M;
end
plot(P_c_array, PR_do_simu_array, 'gs', 'LineWidth',1.0);
hold on;


% case 3
R_d = 1;       % Mbps

% 理论值
PR_do_theo_array = zeros(1,numel(P_c_array));
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
    A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    PR_do_theo_array(idx) = PR_do;
end
plot(P_c_array, PR_do_theo_array, 'b.-', 'LineWidth',1.0);
hold on;

% 仿真值
PR_do_simu_array = zeros(1,numel(P_c_array));
M = 10000;
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    sum_do = 0;     % D2D中断总次数
    for m = 1:M
        % 生成服从均匀分布的随机功率
        P_d = dbm2w(P_dmax) * rand();
        % 生成信道
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % 计算SINR
        SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
        % 计算速率
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
    end
    PR_do_simu_array(idx) = sum_do / M;
end
plot(P_c_array, PR_do_simu_array, 'bs', 'LineWidth',1.0);

grid on;
set(gca,'FontName','Times New Roman');      % 设置坐标轴字体
xlabel('Transmit power of $\mathrm{CUE}_t$, $P_c$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Connection outage probability of D2D link','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend('$R_d$=0.1 Mbps, Theory','$R_d$=0.1 Mbps, Simulation', ...
    '$R_d$=0.5 Mbps, Theory','$R_d$=0.5 Mbps, Simulation', ...
    '$R_d$=1 Mbps, Theory','$R_d$=1 Mbps, Simulation');
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);





