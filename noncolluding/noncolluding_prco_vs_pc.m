%% �����жϸ��� vs P_c
% case 1
clear;
lambda = 1;      % �ŵ�ϵ������
sigma2 = -10;    % dBm����������
P_c_array = 0:2:40;         % dBm
P_dmax = 30;     % dBm
R_c = 0.1;       % Mbps     % case 1

% ����ֵ
PR_co_theo_array = zeros(1,numel(P_c_array));
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
    B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_co_theo_array(idx) = PR_co;
end
plot(P_c_array, PR_co_theo_array, 'r.-', 'LineWidth',1.0);
hold on;

% ����ֵ
PR_co_simu_array = zeros(1,numel(P_c_array));
M = 10000;
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    sum_co = 0;     % �����ж��ܴ���
    for m = 1:M
        % ���ɷ��Ӿ��ȷֲ����������
        P_d = dbm2w(P_dmax) * rand();
        % �����ŵ�
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        % ����SINR
        SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
        % ��������
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
    end
    PR_co_simu_array(idx) = sum_co / M;
end
plot(P_c_array, PR_co_simu_array, 'rs', 'LineWidth',1.0);
hold on;


% case 2
R_c = 0.5;         % Mbps

% ����ֵ
PR_co_theo_array = zeros(1,numel(P_c_array));
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
    B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_co_theo_array(idx) = PR_co;
end
plot(P_c_array, PR_co_theo_array, 'g.-', 'LineWidth',1.0);
hold on;

% ����ֵ
PR_co_simu_array = zeros(1,numel(P_c_array));
M = 10000;
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    sum_co = 0;     % �����ж��ܴ���
    for m = 1:M
        % ���ɷ��Ӿ��ȷֲ����������
        P_d = dbm2w(P_dmax) * rand();
        % �����ŵ�
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        % ����SINR
        SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
        % ��������
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
    end
    PR_co_simu_array(idx) = sum_co / M;
end
plot(P_c_array, PR_co_simu_array, 'gs', 'LineWidth',1.0);
hold on;


% case 3
R_c = 1;         % Mbps

% ����ֵ
PR_co_theo_array = zeros(1,numel(P_c_array));
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
    B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_co_theo_array(idx) = PR_co;
end
plot(P_c_array, PR_co_theo_array, 'b.-', 'LineWidth',1.0);
hold on;

% ����ֵ
PR_co_simu_array = zeros(1,numel(P_c_array));
M = 10000;
for idx = 1:numel(P_c_array)
    P_c = P_c_array(idx);
    sum_co = 0;     % �����ж��ܴ���
    for m = 1:M
        % ���ɷ��Ӿ��ȷֲ����������
        P_d = dbm2w(P_dmax) * rand();
        % �����ŵ�
        h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
        % ����SINR
        SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
        % ��������
        if log2(1+SINR_BS) < R_c
            sum_co = sum_co + 1;
        end
    end
    PR_co_simu_array(idx) = sum_co / M;
end
plot(P_c_array, PR_co_simu_array, 'bs', 'LineWidth',1.0);

grid on;
set(gca,'FontName','Times New Roman');      % ��������������
xlabel('Transmit power of $\mathrm{CUE}_t$, $P_c$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Connection outage probability of cellular link','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend('$R_c$=0.1 Mbps, Theory','$R_c$=0.1 Mbps, Simulation', ...
    '$R_c$=0.5 Mbps, Theory','$R_c$=0.5 Mbps, Simulation', ...
    '$R_c$=1 Mbps, Theory','$R_c$=1 Mbps, Simulation');
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);



