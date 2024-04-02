%% D2D�жϸ��� vs P_dmax
clear;
lambda = 1;      % �ŵ�ϵ������
sigma2 = -10;    % dBm����������

P_c = 15;        % dBm
P_dmax_array = 0:2:40;      % dBm
M = 2*10^4;        % ����ʵ�����

%% case 1
R_d = 0.1;       % Mbps
% ����ֵ
PR_do_theo_array1 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
    A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    PR_do_theo_array1(idx) = PR_do;
end
figure();
plot(P_dmax_array, PR_do_theo_array1, 'b.-', 'LineWidth',1.0);
hold on;

% ����ֵ
PR_do_simu_array1 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_do = 0;     % D2D�ж��ܴ���
    for m = 1:M
        % ���ɷ��Ӿ��ȷֲ����������
        P_d = dbm2w(P_dmax) * rand();
        % �����ŵ�
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % ����SINR
        SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
        % ��������
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
    end
    PR_do_simu_array1(idx) = sum_do / M;
end
plot(P_dmax_array, PR_do_simu_array1, 'bs', 'LineWidth',1.0);
hold on;


%% case 2
R_d = 0.5;       % Mbps
% ����ֵ
PR_do_theo_array2 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
    A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    PR_do_theo_array2(idx) = PR_do;
end
plot(P_dmax_array, PR_do_theo_array2, 'r.--', 'LineWidth',1.0);
hold on;

% ����ֵ
PR_do_simu_array2 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_do = 0;     % D2D�ж��ܴ���
    for m = 1:M
        % ���ɷ��Ӿ��ȷֲ����������
        P_d = dbm2w(P_dmax) * rand();
        % �����ŵ�
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % ����SINR
        SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
        % ��������
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
    end
    PR_do_simu_array2(idx) = sum_do / M;
end
plot(P_dmax_array, PR_do_simu_array2, 'rd', 'LineWidth',1.0);
hold on;


%% case 3
R_d = 1;       % Mbps
% ����ֵ
PR_do_theo_array3 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
    A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    PR_do_theo_array3(idx) = PR_do;
end
plot(P_dmax_array, PR_do_theo_array3, 'g.-.', 'LineWidth',1.0);
hold on;

% ����ֵ
PR_do_simu_array3 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    sum_do = 0;     % D2D�ж��ܴ���
    for m = 1:M
        % ���ɷ��Ӿ��ȷֲ����������
        P_d = dbm2w(P_dmax) * rand();
        % �����ŵ�
        h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
        % ����SINR
        SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
        % ��������
        if log2(1+SINR_DR) < R_d
            sum_do = sum_do + 1;
        end
    end
    PR_do_simu_array3(idx) = sum_do / M;
end
plot(P_dmax_array, PR_do_simu_array3, 'g^', 'LineWidth',1.0);


%%
grid on;
set(gca,'FontName','Times New Roman');      % ��������������
xlabel('Maximum transmit power of DT, $P_d^{\mathrm{max}}$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Connection outage probability of D2D link','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend('$R_d$=0.1 Mbps, Theory','$R_d$=0.1 Mbps, Simulation', ...
    '$R_d$=0.5 Mbps, Theory','$R_d$=0.5 Mbps, Simulation', ...
    '$R_d$=1 Mbps, Theory','$R_d$=1 Mbps, Simulation');
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);





