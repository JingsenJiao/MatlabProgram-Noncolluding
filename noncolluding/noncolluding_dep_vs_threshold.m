%% ��������ߵļ�������� vs �����ֵ
clear;
lambda = 1;      % �ŵ�ϵ������
sigma2 = 10^(-4);      % W����������
P_c = 2;         % W
P_dmax = 15;     % W
threshold_array = 0:0.2:20;
M = 10*10^3;

% h_CTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
% h_DTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
h_CTk = 1;
h_DTk = 1;

% ��������ߵļ�������� vs �����ֵ������ֵ
PFA_theo_array = zeros(1,numel(threshold_array));
PMD_theo_array = zeros(1,numel(threshold_array));
DEP_theo_array = zeros(1,numel(threshold_array));
for idx = 1:numel(threshold_array)
    threshold = threshold_array(idx);
    phi_1 = P_dmax*abs(h_DTk)^2 + sigma2;
    phi_2 = P_c*abs(h_CTk)^2 + sigma2;
    phi_3 = P_dmax*abs(h_DTk)^2 + P_c*abs(h_CTk)^2 + sigma2;
    % FA
    if threshold < sigma2
        PFA = 1;
    elseif threshold <= phi_1
        PFA = 1 - (threshold-sigma2)/(P_dmax*abs(h_DTk)^2);
    else
        PFA = 0;
    end
    PFA_theo_array(idx) = PFA;
    
    % MD
    if threshold < phi_2
        PMD = 0;
    elseif threshold <= phi_3
        PMD = (threshold-phi_2)/(P_dmax*abs(h_DTk)^2);
    else
        PMD = 1;
    end
    PMD_theo_array(idx) = PMD;
    % DEP
    DEP_theo_array(idx) = PFA + PMD;
end


plot(threshold_array, PFA_theo_array, 'g.-', 'LineWidth', 1.0);
hold on;
plot(threshold_array, PMD_theo_array, 'b.-', 'LineWidth', 1.0);
hold on;
plot(threshold_array, DEP_theo_array, 'r.-', 'LineWidth', 1.0);
hold on;

%%
% ��������ߵļ�������� vs �����ֵ������ֵ
PFA_simu_array = zeros(1,numel(threshold_array));
PMD_simu_array = zeros(1,numel(threshold_array));
DEP_simu_array = zeros(1,numel(threshold_array));

for idx = 1:numel(threshold_array)
    threshold = threshold_array(idx);
    
    % FA �龯
    sum_fa = 0;     % �龯���ܴ���
    for m = 1:M
        P_d = P_dmax * rand();   % ���Ӿ��ȷֲ����������
        % ���ʱȽ�
        P_Yk = P_d*abs(h_DTk)^2 + sigma2;
        if P_Yk >= threshold
            sum_fa = sum_fa + 1;
        end
    end
    
    % MD ©��
    sum_md = 0;     % ©����ܴ���
    for m = 1:M
        P_d = P_dmax * rand();   % ���Ӿ��ȷֲ����������
        % ���ʱȽ�
        P_Yk = P_c*abs(h_CTk)^2 + P_d*abs(h_DTk)^2 + sigma2;
        if P_Yk < threshold
            sum_md = sum_md + 1;
        end
    end
    PFA_simu_array(idx) = (sum_fa) / M;
    PMD_simu_array(idx) = (sum_md) / M;
    DEP_simu_array(idx) = (sum_fa + sum_md) / M;
end

plot(threshold_array, PFA_simu_array, 'go', 'LineWidth', 1.0);
hold on;
plot(threshold_array, PMD_simu_array, 'b^', 'LineWidth', 1.0);
hold on;
plot(threshold_array, DEP_simu_array, 'rs', 'LineWidth', 1.0);

grid on;






