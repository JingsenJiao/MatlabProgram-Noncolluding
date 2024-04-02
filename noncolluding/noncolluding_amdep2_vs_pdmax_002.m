%% ��������ߵ�ƽ����С��������� vs P_dmax
% �� �� W ����

clear;
lambda = 1;      % �ŵ�ϵ������
sigma2 = 1;      % W����������
P_c = 5;         % W                  % case 1
P_dmax_array = 1:1:15;     % W


% ƽ����С��������� vs P_dmax������ֵ
AMDEP_theo_array = zeros(1,numel(P_dmax_array));
% AMDEP_theo_array2 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
%     s = (P_c*lambda) / (P_dmax*lambda + P_c*lambda);
%     AMDEP = -s^2 + s*log(s) + 1;
    AMDEP = (P_dmax)/(P_dmax+P_c) + (P_dmax*P_c)/(P_dmax+P_c)^2 + P_c/(P_dmax+P_c)*log(P_c/(P_dmax+P_c));
%     AMDEP2 = s*log(s) - s + 1;
    AMDEP_theo_array(idx) = AMDEP;
%     AMDEP_theo_array2(idx) = AMDEP2;
end
plot(P_dmax_array, AMDEP_theo_array, 'b.-', 'LineWidth', 1.0);
hold on;
% plot(P_dmax_array, AMDEP_theo_array2, 'b.--', 'LineWidth', 1.0);
% hold on;


% ƽ����С��������� vs P_dmax������ֵ
AMDEP_simu_array = zeros(1,numel(P_dmax_array));
M = 5000;
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    % FA �龯
    sum_fa = 0;     % �龯���ܴ���
    for m = 1:M
        P_d = P_dmax * rand();   % ���Ӿ��ȷֲ����������
        h_CTk = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTk = sqrt(lambda/2)*(randn() + 1i*randn());
        % ȷ�����ż����ֵ
        phi_1 = P_dmax * (abs(h_DTk))^2 + sigma2;
        phi_2 = P_c * (abs(h_CTk))^2 + sigma2;
        if phi_1 < phi_2
            threshold = phi_1;
%             threshold = phi_1 + abs(phi_2-phi_1)/2;
        else
            threshold = phi_2;
%             threshold = phi_2 + abs(phi_2-phi_1)/2;
        end
        % ���ʱȽ�
        P_Yk = P_d*(abs(h_DTk))^2 + sigma2;
        if P_Yk >= threshold
            sum_fa = sum_fa + 1;
        end
    end
    
    % MD ©��
    sum_md = 0;     % ©����ܴ���
    for m = 1:M
        P_d = P_dmax * rand();   % ���Ӿ��ȷֲ����������
        h_CTk = sqrt(lambda/2)*(randn() + 1i*randn());
        h_DTk = sqrt(lambda/2)*(randn() + 1i*randn());
        % ȷ�����ż����ֵ
        phi_1 = P_dmax * (abs(h_DTk))^2 + sigma2;
        phi_2 = P_c * (abs(h_CTk))^2 + sigma2;
        if phi_1 < phi_2
            threshold = phi_1;
%             threshold = phi_1 + abs(phi_2-phi_1)/2;
        else
            threshold = phi_2;
%             threshold = phi_2 + abs(phi_2-phi_1)/2;
        end
        % ���ʱȽ�
        P_Yk = P_c*(abs(h_CTk))^2 + P_d*(abs(h_DTk))^2 + sigma2;
        if P_Yk < threshold
            sum_md = sum_md + 1;
        end
    end

    AMDEP_simu_array(idx) = (sum_fa + sum_md) / M;
end

plot(P_dmax_array, AMDEP_simu_array, 'bs', 'LineWidth', 1.0);

grid on;






