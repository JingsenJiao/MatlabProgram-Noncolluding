%% D2D中断概率 vs R_d

clear;
lambda = 1;      % 信道系数方差
sigma2 = -10;    % dBm，噪声功率
P_c = 15;        % dBm
P_dmax = 30;     % dBm
R_d_array = 0:0.1:20;       % Mbps

% 理论值
PR_do_theo_array = zeros(1,numel(R_d_array));
for idx = 1:numel(R_d_array)
    R_d = R_d_array(idx);
    a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
    A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    PR_do_theo_array(idx) = PR_do;
end
figure();
plot(R_d_array, PR_do_theo_array, 'r.-', 'LineWidth',1.0);
hold on;







