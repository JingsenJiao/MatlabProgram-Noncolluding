%% ��������ߵ�ƽ����С��������� vs P_dmax
% �����׵�ͼ3�Աȣ�Covert Communications with A Full-Duplex Receiver over Wireless Fading Channels

clear;
lambda = 1;      % �ŵ�ϵ������
sigma2 = -10;    % dBm����������

P_dmax_array = -20:2:20;       % dBm
M = 10^4;        % ��ʵ�����


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
set(gca,'FontName','Times New Roman');      % ��������������
xlabel('Maximum transmit power of $\mathrm{DUE}_t$, $P_d^{\mathrm{max}}$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Average detection error probability, $\overline{\xi_{k}^*}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend('$P_c$=-20 dBm', '$P_c$=-10 dBm', '$P_c$=0 dBm', '$P_c$=10 dBm', '$P_c$=20 dBm');
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);





