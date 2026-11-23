function [f_WMMSE, C_WMMSE] = WMMSE_ver6(SNR_dB,K,N,H_hat,para)

signal_power = 10^(SNR_dB/10);

f_c = sum(H_hat, 2)/K;

F = zeros(N,K);
for idx_k = 1:K
    F(:, idx_k) = H_hat(:, idx_k);
    F(:, idx_k) = pinv(H_hat*H_hat' + eye*para.sigma^2/signal_power)*H_hat(:,idx_k);
end

curr_rate = 0;
t = 0;
threshold = 0.1;

% f_c = f_c/norm(f_c);
% F = F/norm(F); % normalization 효과 있는지 확인

while 1
    t = t + 1;

    T_c = zeros(1,K);
    T = zeros(1,K);
    I_c = zeros(1,K);
    I = zeros(1,K);
    mmse_c = zeros(1,K);
    mmse = zeros(1,K);
    g_c = zeros(1,K);
    g = zeros(1,K);
    u_c = zeros(1,K);
    u = zeros(1,K);

    %% T, I calculation
    for idx = 1:K
        sum_I = 0;
        for idx_inner = 1:K
            if idx_inner ~= idx
                sum_I = sum_I + abs(H_hat(:,idx)'*F(:,idx_inner))^2;
            end
        end
       
        I_c(idx) = abs(H_hat(:,idx)'*F(:,idx))^2 + sum_I + para.sigma^2;
        T_c(idx) = abs(H_hat(:,idx)'*f_c)^2 + I_c(idx);

        I(idx) = sum_I + para.sigma^2;
        T(idx) = I_c(idx);
    end

    %% mmse & g calculation
    for idx = 1:K
        mmse_c(idx) = I_c(idx)/T_c(idx);
        mmse(idx) = I(idx)/T(idx);

        g_c(idx) = f_c'*H_hat(:,idx)/T_c(idx);
        g(idx) = F(:,idx)'*H_hat(:,idx)/T(idx);

        u_c(idx) = 1/mmse_c(idx);
        u(idx) = 1/mmse(idx);
    end

    cvx_begin quiet

    variable F(N,K) complex
    variable f_c(N) complex
    variable C(K)
    variable R

    maximize R

    subject to
    for idx = 1:K
        sum_I = 0;
        for idx_inner = 1:K
            if idx_inner ~= idx
                sum_I = sum_I + pow_abs(H_hat(:,idx)'*F(:,idx_inner),2);
            end
        end
        
        1 - (u(idx)*(pow_abs(g(idx),2)*(pow_abs(H_hat(:,idx)'*F(:,idx),2) + sum_I + para.sigma^2) - 2*real(g(idx)*H_hat(:,idx)'*F(:,idx))+1) - log(u(idx))/log(2)) + C(idx) >= R
        1 - (u_c(idx)*(pow_abs(g_c(idx),2)*(pow_abs(H_hat(:,idx)'*f_c,2) + pow_abs(H_hat(:,idx)'*F(:,idx),2) + sum_I + para.sigma^2) - 2*real(g_c(idx)*H_hat(:,idx)'*f_c)+1) - log(u_c(idx))/log(2)) >= sum(C)
        C(idx) >= 0
    end

    F_combined = [f_c F];

    norm(F_combined, 'fro') <= sqrt(signal_power)

    cvx_end

    if strcmp(cvx_status, 'Infeasible')
        f_c = sum(H_hat, 2)/K;
        F = zeros(N,K);
        for idx_k = 1:K
            F(:, idx_k) = H_hat(:, idx_k);
        end
        continue
    end

    prev_rate = curr_rate;
    curr_rate = R;

    if abs(prev_rate - curr_rate) < threshold
        break;
    end

    if t > 10
        t = 0;
        % threshold = threshold * 1.5;
    end
end

f_WMMSE = reshape([f_c F],[N*(K+1),1]);
f_WMMSE = f_WMMSE/norm(f_WMMSE,2);
C_WMMSE = C;

end