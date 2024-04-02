%% 平均隐蔽速率 vs P_c
% 按照方法一的检错率作为隐蔽性约束
clear;
lambda = 1;      % 信道系数方差
sigma2 = -10;    % dBm，噪声功率
K = 3;           % 监测者数量
P_c_array = -20:2:40;        % dBm
P_dmax = 10;     % dBm

R_c = 1;       % Mbps
R_d = 1;       % Mbps

M = 10^4;        % 仿真实验次数


% 理论值
actr_theo_array = zeros(1,numel(P_c_array));
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
    B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
    A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    actr_theo_array(idx) = R_c * (1-PR_co) * (1-PR_do);
end

plot(P_c_array, actr_theo_array, 'r.-', 'LineWidth', 1.0);
hold on;


% 仿真值
actr_simu_array = zeros(1,numel(P_c_array));
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    sum_co = 0;     % 蜂窝中断总次数
    sum_do = 0;     % D2D中断总次数
    for m = 1:M
        % 生成服从均匀分布的随机功率
        P_d = dbm2w(P_dmax) * rand();
        % 生成信道
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % 计算SINR
        SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
        SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
        % 比较速率
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
    end
    actr_simu_array(idx) = R_c * (1 - sum_co/M) * (1 - sum_do/M);
end
plot(P_c_array, actr_simu_array, 'rs', 'LineWidth', 1.0);


grid on;
set(gca,'FontName','Times New Roman');      % 设置坐标轴字体
xlabel('Transmit power of $\mathrm{CUE}_t$ $P_c$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Average covert rate (Mbps)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend("$R_c$="+num2str(R_c)+" Mbps, "+"$R_d$="+num2str(R_d)+" Mbps, Theory", ...
                "$R_c$="+num2str(R_c)+" Mbps, "+"$R_d$="+num2str(R_d)+" Mbps, Simulation");
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);






