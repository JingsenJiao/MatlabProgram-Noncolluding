%% 平均最小检测错误概率 vs P_dmax
% case 1
clear;
lambda = 1;      % 信道系数方差
K = 3;           % 监测者总数
sigma2 = -10;    % dBm，噪声功率
P_c = 1;         % dBm                  % case 1
P_dmax_array = 0:2:40;     % dBm
M = 10*10^3;

% 平均最小检测错误概率 vs P_dmax，理论值
AMDEP_theo_array1 = zeros(1,numel(P_dmax_array));
AMDEP_theo_array2 = zeros(1,numel(P_dmax_array));
AMDEP_theo_array3 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    AMDEP1 = (t/(t+1))^K - (K*t^K)/(K+1) * hypergeom([K+1,K+1],K+2,-t);                 % 参考文献
    AMDEP2 = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);   % 自己的
    AMDEP3 = (t/(t+1))^(2*K) - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
    AMDEP_theo_array1(idx) = AMDEP1;
    AMDEP_theo_array2(idx) = AMDEP2;
    AMDEP_theo_array3(idx) = AMDEP3;
end
plot(P_dmax_array, AMDEP_theo_array1, 'b.-', 'LineWidth',1.0);
hold on;
plot(P_dmax_array, AMDEP_theo_array2, 'b.--', 'LineWidth',1.0);
hold on;
plot(P_dmax_array, AMDEP_theo_array3, 'b.-.', 'LineWidth',1.0);
hold on;

% array=[];
% 平均最小检测错误概率 vs P_dmax，仿真值
AMDEP_simu_array = zeros(1,numel(P_dmax_array));  % 仿真值

for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    % FA 虚警
    sum_fa = 0;     % 虚警的总次数
%     threshold_array = zeros(1,M);
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % 服从均匀分布的随机功率
        h_CTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));    % K 个监测者的信道系数
        h_DTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % 最优的监测者
        % 确定最优检测阈值
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        if phi_1 < phi_2
            threshold = phi_1;
        else
            threshold = phi_2;
        end
%         threshold_array(m) = threshold;
        % 功率比较
        P_Yk = P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk >= threshold
            sum_fa = sum_fa + 1;
        end
    end
%     array(numel(array)+1) = mean(threshold_array);
    
    % MD 漏检
    sum_md = 0;     % 漏检的总次数
%     threshold_array = zeros(1,M);
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % 服从均匀分布的随机功率
        h_CTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));    % K 个监测者的信道系数
        h_DTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % 最优的监测者
        % 确定最优检测阈值
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        if phi_1 < phi_2
            threshold = phi_1;
        else
            threshold = phi_2;
        end
%         threshold_array(m) = threshold;
        % 功率比较
        P_Yk = dbm2w(P_c)*abs(h_CTk(k_opt))^2 + P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk < threshold
            sum_md = sum_md + 1;
        end
    end
%     array(numel(array)+1) = mean(threshold_array);

    AMDEP_simu_array(idx) = (sum_fa + sum_md) / M;
end

plot(P_dmax_array, AMDEP_simu_array, 'bs', 'LineWidth',1.0);
hold on;


% case 2
P_c = 10;         % dBm

% 平均最小检测错误概率 vs P_dmax，理论值
AMDEP_theo_array1 = zeros(1,numel(P_dmax_array));
AMDEP_theo_array2 = zeros(1,numel(P_dmax_array));
AMDEP_theo_array3 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    AMDEP1 = (t/(t+1))^K - (K*t^K)/(K+1) * hypergeom([K+1,K+1],K+2,-t);                 % 参考文献
    AMDEP2 = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);   % 自己的
    AMDEP3 = (t/(t+1))^(2*K) - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
    AMDEP_theo_array1(idx) = AMDEP1;
    AMDEP_theo_array2(idx) = AMDEP2;
    AMDEP_theo_array3(idx) = AMDEP3;
end
plot(P_dmax_array, AMDEP_theo_array1, 'r.-', 'LineWidth',1.0);
hold on;
plot(P_dmax_array, AMDEP_theo_array2, 'r.--', 'LineWidth',1.0);
hold on;
plot(P_dmax_array, AMDEP_theo_array3, 'r.-.', 'LineWidth',1.0);
hold on;


% 平均最小检测错误概率 vs P_dmax，仿真值
AMDEP_simu_array = zeros(1,numel(P_dmax_array));  % 仿真值

for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    % FA 虚警
    sum_fa = 0;     % 虚警的总次数
%     threshold_array = zeros(1,M);
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % 服从均匀分布的随机功率
        h_CTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));    % K 个监测者的信道系数
        h_DTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % 最优的监测者
        % 确定最优检测阈值
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        if phi_1 < phi_2
            threshold = phi_1;
        else
            threshold = phi_2;
        end
%         threshold = 0.0058;
%         threshold_array(m) = threshold;
        % 功率比较
        P_Yk = P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk >= threshold
            sum_fa = sum_fa + 1;
        end
    end
%     array(numel(array)+1) = mean(threshold_array);
    
    % MD 漏检
    sum_md = 0;     % 漏检的总次数
%     threshold_array = zeros(1,M);
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % 服从均匀分布的随机功率
        h_CTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));    % K 个监测者的信道系数
        h_DTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % 最优的监测者
        % 确定最优检测阈值
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        if phi_1 < phi_2
            threshold = phi_1;
        else
            threshold = phi_2;
        end
%         threshold = 0.0058;
%         threshold_array(m) = threshold;
        % 功率比较
        P_Yk = dbm2w(P_c)*abs(h_CTk(k_opt))^2 + P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk < threshold
            sum_md = sum_md + 1;
        end
    end
%     array(numel(array)+1) = mean(threshold_array);

    AMDEP_simu_array(idx) = (sum_fa + sum_md) / M;
end

plot(P_dmax_array, AMDEP_simu_array, 'rs', 'LineWidth',1.0);

grid on;
set(gca,'FontName','Times New Roman');      % 设置坐标轴字体
xlabel('Maximum transmit power of $\mathrm{DUE}_t$, $P_d^{\mathrm{max}}$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Average minimum detection error probability, $\overline{\xi_{k^*}^*}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend('$P_c$=5 dBm, Theory ref','$P_c$=5 dBm, Theory','$P_c$=5 dBm, Theory 2','$P_c$=5 dBm, Simulation', ...
    '$P_c$=10 dBm, Theory ref','$P_c$=10 dBm, Theory','$P_c$=10 dBm, Theory 2','$P_c$=10 dBm, Simulation');
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);




