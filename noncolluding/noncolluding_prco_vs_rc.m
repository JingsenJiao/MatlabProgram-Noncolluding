%% 蜂窝中断概率 vs R_c

clear;
lambda = 1;      % 信道系数方差
sigma2 = -10;    % dBm，噪声功率
P_c = 15;        % dBm
P_dmax = 30;     % dBm
R_c_array = 0:0.1:20;       % Mbps

% 理论值
PR_co_theo_array = zeros(1,numel(R_c_array));
for idx = 1:numel(R_c_array)
    R_c = R_c_array(idx);
    b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
    B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_co_theo_array(idx) = PR_co;
end
figure();
plot(R_c_array, PR_co_theo_array, 'r.-', 'LineWidth',1.0);
hold on;





