%% ���ƽ���������� vs ����Լ��
% ���շ���һ�ļ������Ϊ������Լ��
% �ԱȲ�ͬ��K

clear;
lambda = 1;      % �ŵ�ϵ������
sigma2 = -10;    % dBm����������

P_c = 5;         % dBm
omega = 40;      % dBm��pdmax������
P_dmax_array = 0:0.5:omega;     % dBm

R_c = 0.5;       % Mbps
R_d = 0.1;       % Mbps

epsilon_array = 0:0.02:0.3;   % ������Լ��
M = 0*10^3;        % ����ʵ�����

%% case 1
% ����ֵ
K = 1;
max_actr_theo_array = zeros(1,numel(epsilon_array));
for idx = 1:numel(epsilon_array)
    epsilon = epsilon_array(idx);
    actr_theo_array = zeros(1,numel(P_dmax_array));
    for ipd = 1:numel(P_dmax_array)
        P_dmax = P_dmax_array(ipd);
        b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
        B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
        a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
        A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
        PR_co = 1 - exp(-b)*(log(1+B))/(B);
        PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
        t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
        AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
%         AMDEP = (t/(t+1))^K - (K*t^K)/(K+1) * hypergeom([K+1,K+1],K+2,-t);      % �ο�����
        if AMDEP >= (1-epsilon)
            actr_theo_array(ipd) = R_c * (1-PR_co) * (1-PR_do);
        else
            actr_theo_array(ipd) = 0;
        end
    end
    max_actr_theo_array(idx) = max(actr_theo_array);
end

figure();
plot(epsilon_array, max_actr_theo_array, 'b.-', 'LineWidth', 1.0);
hold on;

% ����ֵ
max_actr_simu_array = zeros(1,numel(epsilon_array));
for idx = 1:numel(epsilon_array)
    fprintf("- %d / %d \n", idx, numel(epsilon_array));
    epsilon = epsilon_array(idx);
    actr_simu_array = zeros(1,numel(P_dmax_array));
    for ipd = 1:numel(P_dmax_array)
        if mod(ipd,10)==0
            fprintf("-- %d / %d \n", ipd, numel(P_dmax_array));
        end
        P_dmax = P_dmax_array(ipd);
        sum_co = 0;         % �����ж��ܴ���
        sum_do = 0;         % D2D�ж��ܴ���
        sum_covert = 0;     % ����������Լ�����ܴ���
        for m = 1:M
            % ���ɷ��Ӿ��ȷֲ����������
            P_d = dbm2w(P_dmax) * rand();
            % �����ŵ�
            h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
            h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
            h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
            h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
            % ����SINR
            SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
            SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
            % �Ƚ�����
            if log2(1+SINR_BS) < R_c
                sum_co = sum_co + 1;
            end
            if log2(1+SINR_DR) < R_d
                sum_do = sum_do + 1;
            end
            
            % ����ߵļ��
            sum_fa = 0;
            sum_md = 0;
            for n = 1:M
                P_d_in = dbm2w(P_dmax) * rand();
                h_CTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));    % K ������ߵ��ŵ�ϵ��
                h_DTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));
                mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
                [~,k_opt] = max(mu_k);          % ���ŵļ����
                % ȷ�����ż����ֵ
                phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
                phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
                threshold = min(phi_1, phi_2);
                % FA
                P_Yk_0 = P_d_in*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
                if P_Yk_0 >= threshold
                    sum_fa = sum_fa + 1;
                end
                % MD
                P_Yk_1 = dbm2w(P_c)*abs(h_CTk(k_opt))^2 + P_d_in*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
                if P_Yk_1 < threshold
                    sum_md = sum_md + 1;
                end
            end
            if (sum_fa+sum_md)/M >= 1-epsilon
                sum_covert = sum_covert + 1;
            end
        end     % ������������
        actr_simu_array(ipd) = R_c * (1 - sum_co/M) * (1 - sum_do/M) * sum_covert/M;
    end         % ���� epsilon �ı���
    max_actr_simu_array(idx) = max(actr_simu_array);
end
plot(epsilon_array, max_actr_simu_array, 'bs', 'LineWidth', 1.0);
hold on;


%% case 2
% ����ֵ
K = 3;
max_actr_theo_array = zeros(1,numel(epsilon_array));
for idx = 1:numel(epsilon_array)
    epsilon = epsilon_array(idx);
    actr_theo_array = zeros(1,numel(P_dmax_array));
    for ipd = 1:numel(P_dmax_array)
        P_dmax = P_dmax_array(ipd);
        b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
        B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
        a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
        A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
        PR_co = 1 - exp(-b)*(log(1+B))/(B);
        PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
        t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
        AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
%         AMDEP = (t/(t+1))^K - (K*t^K)/(K+1) * hypergeom([K+1,K+1],K+2,-t);    % �ο�����
        if AMDEP >= (1-epsilon)
            actr_theo_array(ipd) = R_c * (1-PR_co) * (1-PR_do);
        else
            actr_theo_array(ipd) = 0;
        end
    end
    max_actr_theo_array(idx) = max(actr_theo_array);
end

plot(epsilon_array, max_actr_theo_array, 'r.--', 'LineWidth', 1.0);
hold on;

% ����ֵ
max_actr_simu_array = zeros(1,numel(epsilon_array));
for idx = 1:numel(epsilon_array)
    fprintf("- %d / %d \n", idx, numel(epsilon_array));
    epsilon = epsilon_array(idx);
    actr_simu_array = zeros(1,numel(P_dmax_array));
    for ipd = 1:numel(P_dmax_array)
        if mod(ipd,10)==0
            fprintf("-- %d / %d \n", ipd, numel(P_dmax_array));
        end
        P_dmax = P_dmax_array(ipd);
        sum_co = 0;         % �����ж��ܴ���
        sum_do = 0;         % D2D�ж��ܴ���
        sum_covert = 0;     % ����������Լ�����ܴ���
        for m = 1:M
            % ���ɷ��Ӿ��ȷֲ����������
            P_d = dbm2w(P_dmax) * rand();
            % �����ŵ�
            h_CTBS = sqrt(lambda/2)*(randn() + 1i*randn());
            h_DTBS = sqrt(lambda/2)*(randn() + 1i*randn());
            h_CTDR = sqrt(lambda/2)*(randn() + 1i*randn());
            h_DTDR = sqrt(lambda/2)*(randn() + 1i*randn());
            % ����SINR
            SINR_BS = (dbm2w(P_c)*abs(h_CTBS)^2) / (P_d*abs(h_DTBS)^2 + dbm2w(sigma2));
            SINR_DR = (P_d*abs(h_DTDR)^2) / (dbm2w(P_c)*abs(h_CTDR)^2 + dbm2w(sigma2));
            % �Ƚ�����
            if log2(1+SINR_BS) < R_c
                sum_co = sum_co + 1;
            end
            if log2(1+SINR_DR) < R_d
                sum_do = sum_do + 1;
            end
            
            % ����ߵļ��
            sum_fa = 0;
            sum_md = 0;
            for n = 1:M
                P_d_in = dbm2w(P_dmax) * rand();
                h_CTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));    % K ������ߵ��ŵ�ϵ��
                h_DTk = sqrt(lambda/2).*(randn(1,K) + 1i*randn(1,K));
                mu_k = abs(h_CTk).^2 ./ abs(h_DTk).^2;
                [~,k_opt] = max(mu_k);          % ���ŵļ����
                % ȷ�����ż����ֵ
                phi_1 = dbm2w(P_dmax) * abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
                phi_2 = dbm2w(P_c) * abs(h_CTk(k_opt))^2 + dbm2w(sigma2);
                threshold = min(phi_1, phi_2);
                % FA
                P_Yk_0 = P_d_in*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
                if P_Yk_0 >= threshold
                    sum_fa = sum_fa + 1;
                end
                % MD
                P_Yk_1 = dbm2w(P_c)*abs(h_CTk(k_opt))^2 + P_d_in*abs(h_DTk(k_opt))^2 + dbm2w(sigma2);
                if P_Yk_1 < threshold
                    sum_md = sum_md + 1;
                end
            end
            if (sum_fa+sum_md)/M >= 1-epsilon
                sum_covert = sum_covert + 1;
            end
        end     % ������������
        actr_simu_array(ipd) = R_c * (1 - sum_co/M) * (1 - sum_do/M) * sum_covert/M;
    end         % ���� epsilon �ı���
    max_actr_simu_array(idx) = max(actr_simu_array);
end
plot(epsilon_array, max_actr_simu_array, 'rd', 'LineWidth', 1.0);



%%
grid on;
set(gca,'FontName','Times New Roman');      % ��������������
xlabel('Covertness requirement, $\varepsilon$','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Maximum average covert rate (Mbps)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend("K=1, Theory", "K=1, Simulation", ...
                "K=3, Theory", "K=3, Simulation");
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10,'Location','Best');



