%% 蜂窝中断概率 vs P_dmax
clear;
lambda = 1;      % 信道系数方差
sigma2 = -10;    % dBm，噪声功率

P_c = 15;        % dBm
P_dmax_array = 0:2:40;      % dBm
M = 2*10^4;        % 仿真实验次数

%% case 1
R_c = 0.1;       % Mbps
% 理论值
PR_co_theo_array1 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
    B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_co_theo_array1(idx) = PR_co;
end
figure();
plot(P_dmax_array, PR_co_theo_array1, 'b.-', 'LineWidth',1.0);
hold on;

% 仿真值
PR_co_simu_array1 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_co = 0;     % 蜂窝中断总次数
    for m = 1:M
        % 生成服从均匀分布的随机功率
        P_d = dbm2w(P_dmax) * rand();
        % 生成信道
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        % 计算SINR
        SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
        % 计算速率
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
    end
    PR_co_simu_array1(idx) = sum_co / M;
end
plot(P_dmax_array, PR_co_simu_array1, 'bs', 'LineWidth',1.0);
hold on;


%% case 2
R_c = 0.5;       % Mbps
% 理论值
PR_co_theo_array2 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
    B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_co_theo_array2(idx) = PR_co;
end
plot(P_dmax_array, PR_co_theo_array2, 'r.--', 'LineWidth',1.0);
hold on;

% 仿真值
PR_co_simu_array2 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_co = 0;     % 蜂窝中断总次数
    for m = 1:M
        % 生成服从均匀分布的随机功率
        P_d = dbm2w(P_dmax) * rand();
        % 生成信道
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        % 计算SINR
        SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
        % 计算速率
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
    end
    PR_co_simu_array2(idx) = sum_co / M;
end
plot(P_dmax_array, PR_co_simu_array2, 'rd', 'LineWidth',1.0);
hold on;


%% case 3
R_c = 1;       % Mbps
% 理论值
PR_co_theo_array3 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
    B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_co_theo_array3(idx) = PR_co;
end
plot(P_dmax_array, PR_co_theo_array3, 'g.-.', 'LineWidth',1.0);
hold on;

% 仿真值
PR_co_simu_array3 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_co = 0;     % 蜂窝中断总次数
    for m = 1:M
        % 生成服从均匀分布的随机功率
        P_d = dbm2w(P_dmax) * rand();
        % 生成信道
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        % 计算SINR
        SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
        % 计算速率
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
    end
    PR_co_simu_array3(idx) = sum_co / M;
end
plot(P_dmax_array, PR_co_simu_array3, 'g^', 'LineWidth',1.0);


%%
grid on;
set(gca,'FontName','Times New Roman');      % 设置坐标轴字体
xlabel('Maximum transmit power of DT, $P_d^{\mathrm{max}}$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Connection outage probability of cellular link','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend('$R_c$=0.1 Mbps, Theory','$R_c$=0.1 Mbps, Simulation', ...
    '$R_c$=0.5 Mbps, Theory','$R_c$=0.5 Mbps, Simulation', ...
    '$R_c$=1 Mbps, Theory','$R_c$=1 Mbps, Simulation');
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);



