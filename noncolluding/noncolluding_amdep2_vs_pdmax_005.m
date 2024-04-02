% 实数的信道系数
% FA  MD  同时计算

clear;
lambda = 1;
sigma2 = -10;   % dBm
P_c = 5;        % dBm
P_dmax_array = 0:2:30;      % dBm
M = 10^5;


% 平均最小检测错误概率 vs P_dmax，理论值
AMDEP_theo_array = zeros(1,numel(P_dmax_array));
AMDEP_theo_array2 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    s = (dbm2w(P_c)*lambda) / (dbm2w(P_dmax)*lambda + dbm2w(P_c)*lambda);
    AMDEP = -s^2 + s*log(s) + 1;
    AMDEP2 = s*log(s) - s + 1;
    AMDEP_theo_array(idx) = AMDEP;
    AMDEP_theo_array2(idx) = AMDEP2;
end
plot(P_dmax_array, AMDEP_theo_array, 'b.-', 'LineWidth', 1.0);
hold on;
plot(P_dmax_array, AMDEP_theo_array2, 'b.--', 'LineWidth', 1.0);
hold on;


dep_simu_array = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_fa = 0;
    sum_md = 0;
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();
        h_CTk = sqrt(lambda/2) * (randn(1,1)+(-1)*randn(1,1));
        h_DTk = sqrt(lambda/2) * (randn(1,1)+(-1)*randn(1,1));
        phi_1 = dbm2w(P_dmax)*abs(h_DTk) + sigma2;
        phi_2 = dbm2w(P_c)*abs(h_CTk) + sigma2;
        threshold = (phi_1 + phi_2) / 2;
        % FA
        P_Yk_h0 = P_d*abs(h_DTk) + sigma2;
        if P_Yk_h0 >= threshold
            sum_fa = sum_fa + 1;
        end
        % MD
        P_Yk_h1 = dbm2w(P_c)*abs(h_CTk) + P_d*abs(h_DTk) + sigma2;
        if P_Yk_h1 < threshold
            sum_md = sum_md + 1;
        end
    end

    dep_simu_array(idx) = (sum_fa+sum_md) / M;
end

plot(P_dmax_array, dep_simu_array, 'bs', 'LineWidth', 1.0);
grid on;



