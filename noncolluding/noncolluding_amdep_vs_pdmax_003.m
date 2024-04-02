%% ƽ����С��������� vs P_dmax
% ��ĸ��Ź���
% ���Ļ�ͼ����

clear;
lambda = 1;      % �ŵ�ϵ������
K = 3;           % ���������
sigma2 = -10;    % dBm����������

P_dmax_array = 0:1:40;     % dBm
M = 2*10^4;        % ʵ�����
% M = 10^5;        % ʵ�����


%% case 1
P_c = 1;         % dBm

% ƽ����С��������� vs P_dmax������ֵ
AMDEP_theo_array1 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
%     AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) *
%     hypergeom([K+1,K+1],K+2,-t);  % ����ȷ�Ĺ�ʽ
    AMDEP = (t/(t+1))^K - (K*t^K)/(K+1) * hypergeom([K+1,K+1],K+2,-t);
    AMDEP_theo_array1(idx) = AMDEP;
end
figure();
plot(P_dmax_array, AMDEP_theo_array1, 'b.-', 'LineWidth',1.0);
hold on;


% ƽ����С��������� vs P_dmax������ֵ
AMDEP_simu_array1 = zeros(1,numel(P_dmax_array));  % ����ֵ
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    % FA �龯
    sum_fa = 0;     % �龯���ܴ���
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % ���Ӿ��ȷֲ����������
        h_CTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));    % K ������ߵ��ŵ�ϵ��
        h_DTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % ���ŵļ����
        % ȷ�����ż����ֵ
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        if phi_1 < phi_2
            threshold = phi_1;
        else
            threshold = phi_2;
        end
        % ���ʱȽ�
        P_Yk = P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk >= threshold
            sum_fa = sum_fa + 1;
        end
    end
    
    % MD ©��
    sum_md = 0;     % ©����ܴ���
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % ���Ӿ��ȷֲ����������
        h_CTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));    % K ������ߵ��ŵ�ϵ��
        h_DTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % ���ŵļ����
        % ȷ�����ż����ֵ
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        if phi_1 < phi_2
            threshold = phi_1;
        else
            threshold = phi_2;
        end
        % ���ʱȽ�
        P_Yk = dbm2w(P_c)*abs(h_CTk(k_opt))^2 + P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk < threshold
            sum_md = sum_md + 1;
        end
    end

    AMDEP_simu_array1(idx) = (sum_fa + sum_md) / M;
end

plot(P_dmax_array, AMDEP_simu_array1, 'bs', 'LineWidth',1.0);
hold on;


%% case 2
P_c = 5;         % dBm

% ƽ����С��������� vs P_dmax������ֵ
AMDEP_theo_array2 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    % AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
    AMDEP = (t/(t+1))^K - (K*t^K)/(K+1) * hypergeom([K+1,K+1],K+2,-t);
    AMDEP_theo_array2(idx) = AMDEP;
end
plot(P_dmax_array, AMDEP_theo_array2, 'r.--', 'LineWidth',1.0);
hold on;


% ƽ����С��������� vs P_dmax������ֵ
AMDEP_simu_array2 = zeros(1,numel(P_dmax_array));  % ����ֵ
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    % FA �龯
    sum_fa = 0;     % �龯���ܴ���
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % ���Ӿ��ȷֲ����������
        h_CTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));    % K ������ߵ��ŵ�ϵ��
        h_DTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % ���ŵļ����
        % ȷ�����ż����ֵ
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        if phi_1 < phi_2
            threshold = phi_1;
        else
            threshold = phi_2;
        end
        % ���ʱȽ�
        P_Yk = P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk >= threshold
            sum_fa = sum_fa + 1;
        end
    end
    
    % MD ©��
    sum_md = 0;     % ©����ܴ���
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % ���Ӿ��ȷֲ����������
        h_CTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));    % K ������ߵ��ŵ�ϵ��
        h_DTk = sqrt(lambda/2)*(randn(1,K) + 1i*randn(1,K));
        mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
        [~,k_opt] = max(mu_k);          % ���ŵļ����
        % ȷ�����ż����ֵ
        phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
        if phi_1 < phi_2
            threshold = phi_1;
        else
            threshold = phi_2;
        end
        % ���ʱȽ�
        P_Yk = dbm2w(P_c)*abs(h_CTk(k_opt))^2 + P_d*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
        if P_Yk < threshold
            sum_md = sum_md + 1;
        end
    end

    AMDEP_simu_array2(idx) = (sum_fa + sum_md) / M;
end

plot(P_dmax_array, AMDEP_simu_array2, 'rd', 'LineWidth',1.0);


%%
grid on;
set(gca,'FontName','Times New Roman');      % ��������������
xlabel('Maximum transmit power of DT, $P_d^{\mathrm{max}}$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Average minimum detection error probability, $\overline{\xi_{k^*}^*}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend('$P_c$=1 dBm, Theory', '$P_c$=1 dBm, Simulation', ...
                '$P_c$=5 dBm, Theory', '$P_c$=5 dBm, Simulation');
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);



%% ͬʱ�� K=1��K=3
% clear;
% load('D:\MyFiles\˶ʿ\��ҵ�������\B3-ʵ��\ʵ����-����ı�����\С������ͼ\amdep_vs_pdmax_K1-K3_2���.mat');
% figure();
% plot(K1_P_dmax_array, K1_AMDEP_theo_array1, 'g.-.', 'LineWidth',1.0);
% hold on;
% plot(K1_P_dmax_array, K1_AMDEP_simu_array1, 'g^', 'LineWidth',1.0);
% hold on;
% plot(K1_P_dmax_array, K1_AMDEP_theo_array2, 'm.:', 'LineWidth',1.0);
% hold on;
% plot(K1_P_dmax_array, K1_AMDEP_simu_array2, 'mv', 'LineWidth',1.0);
% hold on;
% plot(K1_P_dmax_array, K3_AMDEP_theo_array1, 'b.-', 'LineWidth',1.0);
% hold on;
% plot(K1_P_dmax_array, K3_AMDEP_simu_array1, 'bs', 'LineWidth',1.0);
% hold on;
% plot(K1_P_dmax_array, K3_AMDEP_theo_array2, 'r.--', 'LineWidth',1.0);
% hold on;
% plot(K1_P_dmax_array, K3_AMDEP_simu_array2, 'rd', 'LineWidth',1.0);
% 
% grid on;
% set(gca,'FontName','Times New Roman');      % ��������������
% xlabel('Maximum transmit power of DT, $P_d^{\mathrm{max}}$ (dBm)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
% ylabel('Average minimum detection error probability, $\overline{\xi_{k^*}^*}$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
% handle = legend('$K$=1, $P_c$=1 dBm, Theory', '$K$=1, $P_c$=1 dBm, Simulation', ...
%                 '$K$=1, $P_c$=5 dBm, Theory', '$K$=1, $P_c$=5 dBm, Simulation', ...
%                 '$K$=3, $P_c$=1 dBm, Theory', '$K$=3, $P_c$=1 dBm, Simulation', ...
%                 '$K$=3, $P_c$=5 dBm, Theory', '$K$=3, $P_c$=5 dBm, Simulation');
% set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);





