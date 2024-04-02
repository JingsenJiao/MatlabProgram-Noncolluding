%% 最大平均隐蔽速率 vs 隐蔽约束

clear;
lambda = 1;      % 信道系数方差
sigma2 = -10;    % dBm，噪声功率

P_c = 5;         % dBm
omega = 40;      % dBm，pdmax的上限
P_dmax_array = 0:0.5:omega;     % dBm

R_c = 0.1;       % Mbps
R_d = 0.1;       % Mbps

epsilon_array = 0:0.05:1;   % 隐蔽性约束


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
grid on;
set(gca,'FontName','Times New Roman');      % 设置坐标轴字体
xlabel('Covertness requirement, $\varepsilon$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Maximum average covert rate (Mbps)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
% handle = legend("Theory", "Simulation");
handle = legend("K=1", "K=2", "K=3");
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10,'Location','Best');




