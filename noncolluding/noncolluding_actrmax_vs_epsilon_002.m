%% 最大平均隐蔽速率 vs 隐蔽约束
% 按照方法一的检错率作为隐蔽性约束
% 对比不同的K
% 功率步长越小，曲线越平滑。epsilon步长越大，曲线越平滑。

clear;
lambda = 1;      % 信道系数方差
sigma2 = -10;    % dBm，噪声功率

P_c = 5;         % dBm
omega = 40;      % dBm，pdmax的上限
P_dmax_array = 0:0.5:omega;     % dBm

R_c = 0.1;       % Mbps
R_d = 0.1;       % Mbps

epsilon_array = 0:0.05:1;   % 隐蔽性约束
M = 0*10^3;        % 仿真实验次数

%% case 1
% 理论值
K = 1;
max_actr_theo_array_1 = zeros(1,numel(epsilon_array));
for idx = 1:numel(epsilon_array)
    epsilon = epsilon_array(idx);
    actr_theo_array = zeros(1,numel(P_dmax_array));
    for ipd = 1:numel(P_dmax_array)
        P_dmax = P_dmax_array(ipd);
        b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
        B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
        a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
        A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
        PR_co = 1 - exp(-b)*(log(1+B))/(B);
        PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
        t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
%         AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
        AMDEP = (t/(t+1))^K - (K*t^K)/(K+1) * hypergeom([K+1,K+1],K+2,-t);
        if AMDEP >= (1-epsilon)
            actr_theo_array(ipd) = R_c * (1-PR_co) * (1-PR_do);
        else
            actr_theo_array(ipd) = 0;
        end
    end
    max_actr_theo_array_1(idx) = max(actr_theo_array);
end

figure();
plot(epsilon_array, max_actr_theo_array_1, 'b.-', 'LineWidth', 1.0);
hold on;

%% case 2
% 理论值
K = 2;
max_actr_theo_array_2 = zeros(1,numel(epsilon_array));
for idx = 1:numel(epsilon_array)
    epsilon = epsilon_array(idx);
    actr_theo_array = zeros(1,numel(P_dmax_array));
    for ipd = 1:numel(P_dmax_array)
        P_dmax = P_dmax_array(ipd);
        b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
        B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
        a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
        A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
        PR_co = 1 - exp(-b)*(log(1+B))/(B);
        PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
        t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
%         AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
        AMDEP = (t/(t+1))^K - (K*t^K)/(K+1) * hypergeom([K+1,K+1],K+2,-t);
        if AMDEP >= (1-epsilon)
            actr_theo_array(ipd) = R_c * (1-PR_co) * (1-PR_do);
        else
            actr_theo_array(ipd) = 0;
        end
    end
    max_actr_theo_array_2(idx) = max(actr_theo_array);
end

plot(epsilon_array, max_actr_theo_array_2, 'r.--', 'LineWidth', 1.0);
hold on;

%% case 3
% 理论值
K = 3;
max_actr_theo_array_3 = zeros(1,numel(epsilon_array));
for idx = 1:numel(epsilon_array)
    epsilon = epsilon_array(idx);
    actr_theo_array = zeros(1,numel(P_dmax_array));
    for ipd = 1:numel(P_dmax_array)
        P_dmax = P_dmax_array(ipd);
        b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
        B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
        a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
        A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
        PR_co = 1 - exp(-b)*(log(1+B))/(B);
        PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
        t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
%         AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
        AMDEP = (t/(t+1))^K - (K*t^K)/(K+1) * hypergeom([K+1,K+1],K+2,-t);
        if AMDEP >= (1-epsilon)
            actr_theo_array(ipd) = R_c * (1-PR_co) * (1-PR_do);
        else
            actr_theo_array(ipd) = 0;
        end
    end
    max_actr_theo_array_3(idx) = max(actr_theo_array);
end

plot(epsilon_array, max_actr_theo_array_3, 'g.-.', 'LineWidth', 1.0);
hold on;


%%
% % 仿真值
% max_actr_simu_array = zeros(1,numel(epsilon_array));
% for idx = 1:numel(epsilon_array)
%     fprintf("- %d / %d \n", idx, numel(epsilon_array));
%     epsilon = epsilon_array(idx);
%     actr_simu_array = zeros(1,numel(P_dmax_array));
%     for ipd = 1:numel(P_dmax_array)
%         fprintf("-- %d / %d \n", ipd, numel(P_dmax_array));
%         P_dmax = P_dmax_array(ipd);
%         sum_co = 0;         % 蜂窝中断总次数
%         sum_do = 0;         % D2D中断总次数
%         sum_covert = 0;     % 满足隐蔽性约束的总次数
%         for m = 1:M
%             % 生成服从均匀分布的随机功率
%             P_d = dbm2w(P_dmax) * rand();
%             % 生成信道
%             h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
%             h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
%             h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
%             h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
%             % 计算SINR
%             SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
%             SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
%             % 比较速率
%             if log2(1+SINR_BS) < R_c
%                 sum_co = sum_co + 1;
%             end
%             if log2(1+SINR_DR) < R_d
%                 sum_do = sum_do + 1;
%             end
%             
%             % 监测者的检测
%             sum_fa = 0;
%             sum_md = 0;
%             for n = 1:M
%                 P_d_in = dbm2w(P_dmax) * rand();
%                 h_CTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));    % K 个监测者的信道系数
%                 h_DTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));
%                 mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
%                 [~,k_opt] = max(mu_k);          % 最优的监测者
%                 % 确定最优检测阈值
%                 phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
%                 phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
%                 threshold = min(phi_1, phi_2);
%                 % FA
%                 P_Yk_0 = P_d_in*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
%                 if P_Yk_0 >= threshold
%                     sum_fa = sum_fa + 1;
%                 end
%                 % MD
%                 P_Yk_1 = dbm2w(P_c)*abs(h_CTk(k_opt))^2 + P_d_in*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
%                 if P_Yk_1 < threshold
%                     sum_md = sum_md + 1;
%                 end
%             end
%             if (sum_fa+sum_md)/M >= 1-epsilon
%                 sum_covert = sum_covert + 1;
%             end
%         end     % 结束最外层仿真
%         actr_simu_array(ipd) = R_c * (1 - sum_co/M) * (1 - sum_do/M) * sum_covert/M;
%     end         
%     max_actr_simu_array(idx) = max(actr_simu_array);
% end           % 结束 epsilon 的遍历
% plot(epsilon_array, max_actr_simu_array, 'bs', 'LineWidth', 1.0);


%%
grid on;
set(gca,'FontName','Times New Roman');      % 设置坐标轴字体
xlabel('Covertness requirement, $\varepsilon$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Maximum average covert rate (Mbps)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
% handle = legend("Theory", "Simulation");
handle = legend("K=1", "K=2", "K=3");
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10,'Location','Best');




