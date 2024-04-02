% 比较方法1、方法2

clear;
lambda = 1;      % 信道系数方差
sigma2 = -10;    % dBm，噪声功率
P_c = 10;         % dBm
P_dmax_array = -20:2:40;     % dBm

K = 2;           % 监测者总数
M = 10^4;        % 仿真实验次数

% 理论值
AMDEP_theo_array1 = zeros(1,numel(P_dmax_array));
AMDEP_theo_array2 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);

    t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    AMDEP1 = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);

    s = (dbm2w(P_c)*lambda) / (dbm2w(P_dmax)*lambda + dbm2w(P_c)*lambda);
    AMDEP2 = -s^2 + s*log(s) + 1;

    AMDEP_theo_array1(idx) = AMDEP1;
    AMDEP_theo_array2(idx) = AMDEP2;
end

plot(P_dmax_array, AMDEP_theo_array1, 'r.-', 'LineWidth', 1.0);
hold on;
plot(P_dmax_array, AMDEP_theo_array2, 'b.--', 'LineWidth', 1.0);
hold on;


% 仿真值，1
AMDEP_simu_array1 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_fa = 0;
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % 服从均匀分布的随机功率
        h_CTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));    % K 个监测者的信道系数
        h_DTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % 最优的监测者
        % 确定最优检测阈值
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        threshold = min(phi_1, phi_2);
        % 功率比较
        P_Yk = P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk >= threshold
            sum_fa = sum_fa + 1;
        end
    end
    
    sum_md = 0;
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % 服从均匀分布的随机功率
        h_CTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));    % K 个监测者的信道系数
        h_DTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % 最优的监测者
        % 确定最优检测阈值
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        threshold = min(phi_1, phi_2);
        % 功率比较
        P_Yk = dbm2w(P_c)*abs(h_CTk(k_opt))^2 + P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk < threshold
            sum_md = sum_md + 1;
        end
    end

    AMDEP_simu_array1(idx) = (sum_fa + sum_md) / M;
end

plot(P_dmax_array, AMDEP_simu_array1, 'ro', 'LineWidth', 1.0);
hold on;


% 仿真值，2
AMDEP_simu_array1 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_fa = 0;
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % 服从均匀分布的随机功率
        h_CTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
        h_DTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
        % 确定最优检测阈值
        phi_1 = dbm2w(P_dmax) * abs(h_DTk)^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk)^2 + dbm2w(sigma2);
        threshold = min(phi_1, phi_2);
        % 功率比较
        P_Yk = P_d*abs(h_DTk)^2 + dbm2w(sigma2);
        if P_Yk >= threshold
            sum_fa = sum_fa + 1;
        end
    end

    sum_md = 0;
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % 服从均匀分布的随机功率
        h_CTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
        h_DTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
        % 确定最优检测阈值
        phi_1 = dbm2w(P_dmax) * abs(h_DTk)^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk)^2 + dbm2w(sigma2);
        threshold = min(phi_1, phi_2);
        % 功率比较
        P_Yk = dbm2w(P_c)*abs(h_CTk)^2 + P_d*abs(h_DTk)^2 + dbm2w(sigma2);
        if P_Yk < threshold
            sum_md = sum_md + 1;
        end
    end

    AMDEP_simu_array1(idx) = (sum_fa + sum_md) / M;
end

plot(P_dmax_array, AMDEP_simu_array1, 'bs', 'LineWidth', 1.0);

grid on;

set(gca,'FontName','Times New Roman');      % 设置坐标轴字体
xlabel('Maximum transmit power of $\mathrm{DUE}_t$ $P_d^{\mathrm{max}}$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Average minimum detection error probability $\overline{\xi_{k}^*}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend('1 Theory','2 Theory','1 Simulation','2 Simulation');
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);



