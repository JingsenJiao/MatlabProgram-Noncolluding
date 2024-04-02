%% 平均隐蔽速率 vs P_dmax
% 按照方法一的检错率作为隐蔽性约束
% 当R_d在7 Mbps以内时，R_c、R_d的取值不影响最优的P_dmax的选择。

clear;
lambda = 1;      % 信道系数方差
sigma2 = 10^(-4);    % W，噪声功率
K = 3;           % 监测者数量

P_dmax_array = 0:0.5:10;     % W

R_c = 0.5;       % Mbps
R_d = 0.1;       % Mbps

epsilon = 0.1;   % 隐蔽性约束
M = 0.9*10^3;        % 仿真实验次数


%% case 1
P_c = 0.01;         % W
% 理论值
actr_theo_array = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    b = ((2^R_c-1)*sigma2) / (lambda*P_c);
    B = ((2^R_c-1)*lambda*P_dmax) / (lambda*P_c);
    a = sigma2 / (lambda*P_c);
    A = ((2^R_d-1)*lambda*P_c) / (lambda*P_dmax);
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    t = (lambda*P_dmax) / (lambda*P_c);
    AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
    if AMDEP >= (1-epsilon)
        actr_theo_array(idx) = R_c * (1-PR_co) * (1-PR_do);
    else
        actr_theo_array(idx) = 0;
    end
end

figure();
plot(P_dmax_array, actr_theo_array, 'b.-', 'LineWidth', 1.0);
hold on;

% 仿真值
actr_simu_array = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    fprintf("%d / %d \n", idx, numel(P_dmax_array));
    P_dmax = P_dmax_array(idx);
    sum_co = 0;         % 蜂窝中断总次数
    sum_do = 0;         % D2D中断总次数
    sum_covert = 0;     % 满足隐蔽性约束的总次数
    for m = 1:M
        if mod(m,1000)==0
            fprintf("%d / %d \n", m, M);
        end
        % 生成服从均匀分布的随机功率
        P_d = P_dmax * rand();
        % 生成信道
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % 计算SINR
        SINR_BS = (P_c*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + sigma2);
        SINR_DR = (P_d*abs(h_DTDR)^2) / (P_c*abs(h_CTDR)^2 + sigma2);
        % 比较速率
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
        
        % 监测者的检测
        sum_fa = 0;
        sum_md = 0;
        for n = 1:M
            P_d_in = P_dmax * rand();
            h_CTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));    % K 个监测者的信道系数
            h_DTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));
            mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
            [~,k_opt] = max(mu_k);          % 最优的监测者
            % 确定最优检测阈值
            phi_1 = P_dmax * abs(h_DTk(k_opt))^2 + sigma2;
            phi_2 = P_c * abs(h_CTk(k_opt))^2 + sigma2;
            threshold = min(phi_1, phi_2);
            % FA
            P_Yk_0 = P_d_in*abs(h_DTk(k_opt))^2 + sigma2;
            if P_Yk_0 >= threshold
                sum_fa = sum_fa + 1;
            end
            % MD
            P_Yk_1 = P_c*abs(h_CTk(k_opt))^2 + P_d_in*abs(h_DTk(k_opt))^2 + sigma2;
            if P_Yk_1 < threshold
                sum_md = sum_md + 1;
            end
        end
        if (sum_fa+sum_md)/M >= 1-epsilon
            sum_covert = sum_covert + 1;
        end

    end
    actr_simu_array(idx) = R_c * (1 - sum_co/M) * (1 - sum_do/M) * sum_covert/M;
end
plot(P_dmax_array, actr_simu_array, 'bs', 'LineWidth', 1.0);
hold on;


%% case 2
P_c = 0.03;         % W
% 理论值
actr_theo_array = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    b = ((2^R_c-1)*sigma2) / (lambda*P_c);
    B = ((2^R_c-1)*lambda*P_dmax) / (lambda*P_c);
    a = sigma2 / (lambda*P_c);
    A = ((2^R_d-1)*lambda*P_c) / (lambda*P_dmax);
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    t = (lambda*P_dmax) / (lambda*P_c);
    AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
    if AMDEP >= (1-epsilon)
        actr_theo_array(idx) = R_c * (1-PR_co) * (1-PR_do);
    else
        actr_theo_array(idx) = 0;
    end
end

plot(P_dmax_array, actr_theo_array, 'r.--', 'LineWidth', 1.0);
hold on;

% 仿真值
actr_simu_array = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    fprintf("%d / %d \n", idx, numel(P_dmax_array));
    P_dmax = P_dmax_array(idx);
    sum_co = 0;         % 蜂窝中断总次数
    sum_do = 0;         % D2D中断总次数
    sum_covert = 0;     % 满足隐蔽性约束的总次数
    for m = 1:M
        if mod(m,1000)==0
            fprintf("%d / %d \n", m, M);
        end
        % 生成服从均匀分布的随机功率
        P_d = P_dmax * rand();
        % 生成信道
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % 计算SINR
        SINR_BS = (P_c*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + sigma2);
        SINR_DR = (P_d*abs(h_DTDR)^2) / (P_c*abs(h_CTDR)^2 + sigma2);
        % 比较速率
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
        
        % 监测者的检测
        sum_fa = 0;
        sum_md = 0;
        for n = 1:M
            P_d_in = P_dmax * rand();
            h_CTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));    % K 个监测者的信道系数
            h_DTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));
            mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
            [~,k_opt] = max(mu_k);          % 最优的监测者
            % 确定最优检测阈值
            phi_1 = P_dmax * abs(h_DTk(k_opt))^2 + sigma2;
            phi_2 = P_c * abs(h_CTk(k_opt))^2 + sigma2;
            threshold = min(phi_1, phi_2);
            % FA
            P_Yk_0 = P_d_in*abs(h_DTk(k_opt))^2 + sigma2;
            if P_Yk_0 >= threshold
                sum_fa = sum_fa + 1;
            end
            % MD
            P_Yk_1 = P_c*abs(h_CTk(k_opt))^2 + P_d_in*abs(h_DTk(k_opt))^2 + sigma2;
            if P_Yk_1 < threshold
                sum_md = sum_md + 1;
            end
        end
        if (sum_fa+sum_md)/M >= 1-epsilon
            sum_covert = sum_covert + 1;
        end

    end
    actr_simu_array(idx) = R_c * (1 - sum_co/M) * (1 - sum_do/M) * sum_covert/M;
end
plot(P_dmax_array, actr_simu_array, 'rd', 'LineWidth', 1.0);



grid on;
set(gca,'FontName','Times New Roman');      % 设置坐标轴字体
xlabel('Maximum transmit power of $\mathrm{DUE}_t$, $P_d^{\mathrm{max}}$ (W)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Average covert rate (Mbps)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
% handle = legend("$R_c$="+num2str(R_c)+" Mbps, "+"$R_d$="+num2str(R_d)+" Mbps, $P_c$=5 dBm, Theory", ...
%                 "$R_c$="+num2str(R_c)+" Mbps, "+"$R_d$="+num2str(R_d)+" Mbps, $P_c$=5 dBm, Simulation", ...
%                 "$R_c$="+num2str(R_c)+" Mbps, "+"$R_d$="+num2str(R_d)+" Mbps, $P_c$=10 dBm, Theory", ...
%                 "$R_c$="+num2str(R_c)+" Mbps, "+"$R_d$="+num2str(R_d)+" Mbps, $P_c$=10 dBm, Simulation");
handle = legend("$P_c$=0.01 W, Theory", ...
                "$P_c$=0.01 W, Simulation", ...
                "$P_c$=0.03 W, Theory", ...
                "$P_c$=0.03 W, Simulation");
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10,'Location','Best');






