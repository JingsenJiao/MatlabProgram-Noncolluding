%% ƽ���������� vs P_dmax
% �� W Ϊ���ʵ�λ
% ���շ���һ�ļ������Ϊ������Լ��
clear;
lambda = 1;      % �ŵ�ϵ������

sigma2 = 10^(-4);    % W����������
K = 3;           % ���������
P_c = 0.1;         % W
P_dmax_array = 0:0.05:1;     % W

R_c = 0.5;       % Mbps
R_d = 1;       % Mbps

M = 8*10^3;        % ����ʵ�����

% ����ֵ
actr_theo_array = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    b = ((2^R_c-1)*sigma2) / (lambda*P_c);
    B = ((2^R_c-1)*lambda*P_dmax) / (lambda*P_c);
    a = sigma2 / (lambda*P_c);
    A = ((2^R_d-1)*lambda*P_c) / (lambda*P_dmax);
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    actr_theo_array(idx) = R_c * (1-PR_co) * (1-PR_do);
end

plot(P_dmax_array, actr_theo_array, 'r.-', 'LineWidth', 1.0);
hold on;


% ����ֵ
actr_simu_array = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_co = 0;         % �����ж��ܴ���
    sum_do = 0;         % D2D�ж��ܴ���
    for m = 1:M
        % ���ɷ��Ӿ��ȷֲ����������
        P_d = P_dmax * rand();
        % �����ŵ�
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % ����SINR
        SINR_BS = (P_c*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + sigma2);
        SINR_DR = (P_d*abs(h_DTDR)^2) / (P_c*abs(h_CTDR)^2 + sigma2);
        % �Ƚ�����
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
    end
    actr_simu_array(idx) = R_c * (1 - sum_co/M) * (1 - sum_do/M);
end
plot(P_dmax_array, actr_simu_array, 'rs', 'LineWidth', 1.0);


grid on;
set(gca,'FontName','Times New Roman');      % ��������������
xlabel('Maximum transmit power of $\mathrm{DUE}_t$ $P_d^{\mathrm{max}}$ (W)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Average covert rate (Mbps)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend("$R_c$="+num2str(R_c)+" Mbps, "+"$R_d$="+num2str(R_d)+" Mbps, Theory", ...
                "$R_c$="+num2str(R_c)+" Mbps, "+"$R_d$="+num2str(R_d)+" Mbps, Simulation");
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);






