%% ƽ���������� vs R_c
% ���շ���һ�ļ������Ϊ������Լ��
clear;
lambda = 1;      % �ŵ�ϵ������
sigma2 = -10;    % dBm����������
K = 3;           % ���������
P_c = 5;        % dBm
P_dmax = 30;     % dBm

R_c_array = 0:0.2:3;       % Mbps
R_d = 0.1;       % Mbps

epsilon = 0.1;   % ������Լ��
M = 20*10^3;        % ����ʵ�������2w����Ȼ�Ƚ���ɢ


% ����ֵ
actr_theo_array = zeros(1,numel(R_c_array));
for idx = 1:numel(R_c_array)
    R_c = R_c_array(idx);
    b = ((2^R_c-1)*dbm2w(sigma2)) / (lambda*dbm2w(P_c));
    B = ((2^R_c-1)*lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    a = dbm2w(sigma2) / (lambda*dbm2w(P_c));
    A = ((2^R_d-1)*lambda*dbm2w(P_c)) / (lambda*dbm2w(P_dmax));
    PR_co = 1 - exp(-b)*(log(1+B))/(B);
    PR_do = 1 - exp(-a*A) + A*(a+1)*expint(a*A) - A*exp(a)*expint(a*(A+1));
    t = (lambda*dbm2w(P_dmax)) / (lambda*dbm2w(P_c));
    AMDEP = (t/(t+1))^K - (K*t^(2*K))/((K+1)*(t+1)^K) * hypergeom([K+1,K+1],K+2,-t);
    if AMDEP >= (1-epsilon)
        actr_theo_array(idx) = R_c * (1-PR_co) * (1-PR_do);
    else
        actr_theo_array(idx) = 0;
    end
end
figure();
plot(R_c_array, actr_theo_array, 'r.-', 'LineWidth', 1.0);
hold on;


% ����ֵ
actr_simu_array = zeros(1,numel(R_c_array));

for idx = 1:numel(R_c_array)
    fprintf("- %d / %d\n", idx, numel(R_c_array));
    R_c = R_c_array(idx);
    sum_co = 0;     % �����ж��ܴ���
    sum_do = 0;     % D2D�ж��ܴ���
    sum_covert = 0;     % ����������Լ�����ܴ���
    for m = 1:M
        if mod(m,1000)==0
            fprintf("-- %d / %d\n", m, M);
        end
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

    end
    actr_simu_array(idx) = R_c * (1 - sum_co/M) * (1 - sum_do/M) * sum_covert/M;
end
plot(R_c_array, actr_simu_array, 'rs', 'LineWidth', 1.0);


grid on;
set(gca,'FontName','Times New Roman');      % ��������������
xlabel('Predetermined transmission rate of cellular link, $R_c$ (Mbps)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
ylabel('Average covert rate (Mbps)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
handle = legend("$R_d$="+num2str(R_d)+" Mbps, Theory", ...
                "$R_d$="+num2str(R_d)+" Mbps, Simulation");
set(handle,'Interpreter','latex','FontName','Times New Roman','FontSize',10);






