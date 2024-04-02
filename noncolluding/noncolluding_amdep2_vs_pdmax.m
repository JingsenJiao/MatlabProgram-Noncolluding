%% ��������ߵ�ƽ����С��������� vs P_dmax
% �����Ź��������ι��ʣ���WΪ��λ����16�����ϡ���0.01������ʱ������������Ǻ�

% digits(32)
clear;
lambda = 1;      % �ŵ�ϵ������
sigma2 = -10;    % dBm����������
P_c = 5;         % dBm
P_dmax_array = -20:2:40;     % dBm
M = 10^4;        % ��ʵ�����


% ƽ����С��������� vs P_dmax������ֵ
AMDEP_theo_array = zeros(1,numel(P_dmax_array));
AMDEP_theo_array2 = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    s = (dbm2w(P_c)*lambda) / (dbm2w(P_dmax)*lambda + dbm2w(P_c)*lambda);
    AMDEP = -s^2 + s*log(s) + 1;
    AMDEP2 = s*log(s) - s + 1;
    AMDEP_theo_array(idx) = AMDEP;
    AMDEP_theo_array2(idx) = AMDEP2;
end
plot(P_dmax_array, AMDEP_theo_array, 'b.-', 'LineWidth', 1.0);
hold on;
plot(P_dmax_array, AMDEP_theo_array2, 'b.--', 'LineWidth', 1.0);
hold on;


% ƽ����С��������� vs P_dmax������ֵ
AMDEP_simu_array = zeros(1,numel(P_dmax_array));
for idx = 1:numel(P_dmax_array)
    P_dmax = P_dmax_array(idx);
    % FA �龯
    sum_fa = 0;     % �龯���ܴ���
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % ���Ӿ��ȷֲ����������
        h_CTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
        h_DTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
        % ȷ�����ż����ֵ
        phi_1 = dbm2w(P_dmax) * abs(h_DTk)^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk)^2 + dbm2w(sigma2);
        threshold = min(phi_1, phi_2);
        % ���ʱȽ�
        P_Yk = P_d*abs(h_DTk)^2 + dbm2w(sigma2);
        if P_Yk >= threshold
            sum_fa = sum_fa + 1;
        end
    end
    
    % MD ©��
    sum_md = 0;     % ©����ܴ���
    for m = 1:M
        P_d = dbm2w(P_dmax) * rand();   % ���Ӿ��ȷֲ����������
        h_CTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
        h_DTk = sqrt(lambda/2)*(randn(1,1) + 1i*randn(1,1));
        % ȷ�����ż����ֵ
        phi_1 = dbm2w(P_dmax) * abs(h_DTk)^2 + dbm2w(sigma2);
        phi_2 = dbm2w(P_c) * abs(h_CTk)^2 + dbm2w(sigma2);
        threshold = min(phi_1, phi_2);
        % ���ʱȽ�
        P_Yk = dbm2w(P_c)*abs(h_CTk)^2 + P_d*abs(h_DTk)^2 + dbm2w(sigma2);
        if P_Yk < threshold
            sum_md = sum_md + 1;
        end
    end

    AMDEP_simu_array(idx) = (sum_fa + sum_md) / M;
end

plot(P_dmax_array, AMDEP_simu_array, 'bs', 'LineWidth', 1.0);

grid on;






