%% 单个监测者的平均最小检测错误概率 vs P_dmax
% 与文献的图3对比：Covert Communications with A Full-Duplex Receiver over Wireless Fading Channels

clear;
lambda = 1;      % 信道系数方差
sigma2 = -10;    % dBm，噪声功率

P_dmax_array = -20:2:20;       % dBm
M = 10^4;        % 总实验次数


P_c_array = -20:10:20;       % dBm
color_array = ["b.-", "m.-", "y.-", "c.-", "g.-"];
for ic = 1:numel(P_c_array)
    P_c = P_c_array(ic);
    AMDEP_theo_array = zeros(1,numel(P_dmax_array));
    for idx = 1:numel(P_dmax_array)
        P_dmax = P_dmax_array(idx);
        s = (dbm2w(P_c)*lambda) / (dbm2w(P_dmax)*lambda + dbm2w(P_c)*lambda);
        AMDEP = -s^2 + s*log(s) + 1;
        AMDEP_theo_array(idx) = AMDEP;
    end
    plot(P_dmax_array, AMDEP_theo_array, color_array(ic), 'LineWidth', 1.0);
    hold on;
end

grid on;
set(gca,'FontName','Times New Roman');      % 设置坐标轴字体
xlabel('Maximum transmit power of $\mathrm{DUE}_t$, $P_d^{\mathrm{max}}$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Average detection error probability, $\overline{\xi_{k}^*}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend('$P_c$=-20 dBm', '$P_c$=-10 dBm', '$P_c$=0 dBm', '$P_c$=10 dBm', '$P_c$=20 dBm');
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);





