function [f_GPI_RS, C_portion, max_idx] = MMF_GPI_RS(lag_multi_range,K,N,H_hat,H_hat_err_cov,SNR_N, para)

cnt_lag = 0;

% matrix computation
for idx_k = 1:K
    h = H_hat(:, idx_k);
    phi = H_hat_err_cov(:, :, idx_k);
    % phi = 0;
    
    A{idx_k}(:,:) = zeros(N*(K+1), N*(K+1));
    B{idx_k}(:,:) = zeros(N*(K+1), N*(K+1));

    for idx_k_inner = 1:1+K
        if idx_k_inner == 1
            A_index = (idx_k_inner-1)*N + 1: (idx_k_inner-1)*N + N;
            A{idx_k}(A_index, A_index) = zeros(N,N);
        else
            A_index = (idx_k_inner-1)*N + 1: (idx_k_inner-1)*N + N;
            A{idx_k}(A_index, A_index) = h*h' + phi;
        end
    end
    
    A{idx_k} = A{idx_k} + eye(N*(K+1), N*(K+1))*(1/SNR_N);
    B{idx_k} = A{idx_k};
    
    B_index = idx_k*N + 1: idx_k*N + N;
    B{idx_k}(B_index, B_index) = eye(N, N)*(1/SNR_N) + phi; % phi 여기 추가됨


    C{idx_k}(:,:) = zeros(N*(K+1), N*(K+1));
    D{idx_k}(:,:) = zeros(N*(K+1), N*(K+1));

    for idx_k_inner = 1:1+K
        C_index = (idx_k_inner-1)*N + 1: (idx_k_inner-1)*N + N;
        C{idx_k}(C_index, C_index) = h*h' + phi;
    end

    C{idx_k} = C{idx_k} + eye(N*(K+1), N*(K+1))*(1/SNR_N);
    D{idx_k} = C{idx_k};

    D_index = 1:N;
    D{idx_k}(D_index, D_index) = eye(N, N)*(1/SNR_N) + phi; % 여기 추가됨

end


%% ZF version
F_fd_MSE = pinv(H_hat*H_hat' + (1/SNR_N)*eye(N, N))*H_hat;

F_RS = [1/K*sum(F_fd_MSE,2), F_fd_MSE];
f_curr = reshape(F_RS,N*(K+1),1);


% %% MRT version
% f_mtx = zeros(N,K);
% 
% for idx_k = 1:K
%     f_mtx(:, idx_k) = H_hat(:, idx_k);
% end
% F_RS = [1/K*sum(f_mtx,2), f_mtx];
% f_curr = reshape(F_RS,N*(K+1),1);

%%
f_curr = f_curr/norm(f_curr);

for lag_multi = lag_multi_range
    cnt_lag = cnt_lag + 1;
    
    t = 0;
    alpha = para.alpha;

    while 1
        t = t + 1;
        f_prev = f_curr;
        f = f_curr;

        for idx_k = 1:K
            rate_private(idx_k) = real(log2(f'*A{idx_k}*f/(f'*B{idx_k}*f)));
            rate_common(idx_k) = real(log2(f'*C{idx_k}*f/(f'*D{idx_k}*f)));
        end

        % first time, we use uniformly distributed C_portion
        % if cnt_lag == 1
        C_portion = min(rate_common)/K * ones(1,K);
        % end

        for idx_k = 1:K
            w_1(idx_k) = exp(-(rate_common(idx_k))/alpha);
            w_2(idx_k) = exp(-(rate_private(idx_k)+C_portion(idx_k))/alpha);
        end

        w_1 = w_1/sum(w_1);
        w_2 = w_2/sum(w_2);

        for idx_k = 1:K
            % common message matrix
            X_kkt_common{idx_k} = lag_multi*w_1(idx_k)*C{idx_k}/(f'*C{idx_k}*f);
            Y_kkt_common{idx_k} = lag_multi*w_1(idx_k)*D{idx_k}/(f'*D{idx_k}*f);
    
            % private + C_portion matrix
            X_kkt_private{idx_k} = w_2(idx_k)*A{idx_k}/(f'*A{idx_k}*f);
            Y_kkt_private{idx_k} = w_2(idx_k)*B{idx_k}/(f'*B{idx_k}*f);
        end
    
        X_kkt = 0;
        Y_kkt = 0;
    
        for idx_k = 1:K
            X_kkt = X_kkt + X_kkt_common{idx_k} + X_kkt_private{idx_k};
            Y_kkt = Y_kkt + Y_kkt_common{idx_k} + Y_kkt_private{idx_k};
        end
    
        for idx_k = 1:K+1
            Y_kkt_part{idx_k} = Y_kkt((idx_k-1)*N + 1:(idx_k)*N,(idx_k-1)*N + 1:(idx_k)*N);
            inv_Y_kkt_part{idx_k} = pinv(Y_kkt_part{idx_k});
        end
    
        inv_Y_kkt = blkdiag(inv_Y_kkt_part{:});
    
        f_curr = inv_Y_kkt*X_kkt*f;
        % f_curr = pinv(Y_kkt)*X_kkt*f;
    
        f_curr = f_curr/norm(f_curr, 2);

        if norm(f_prev-f_curr, 2) < para.th
            break;  
        end

        if t > para.loop_terminate
            alpha = alpha*(1+para.alpha_increase);
    
            t = 0;
            f = randn(N*(K+1),1) + 1i*randn(N*(K+1),1);
            f = f/norm(f,2);
    
            f_curr = f;
        end

    end

    % given f_curr, step 2 is initiated
    f_curr;

    for idx_k = 1:K
        rate_private(idx_k) = real(log2(f_curr'*A{idx_k}*f_curr/(f_curr'*B{idx_k}*f_curr)));
        rate_common(idx_k) = real(log2(f_curr'*C{idx_k}*f_curr/(f_curr'*D{idx_k}*f_curr)));
    end

    %% Stage2 (f_curr is given and fixed)
    
    [sorted_R, idx_R] = sort(rate_private);
    reverse_idx(idx_R) = 1:K;

    L_indicator = zeros(1, para.K);

    for idx_i = 1:K-1
        for idx_j = 1:idx_i
            L_indicator(idx_i) = L_indicator(idx_i) + sorted_R(idx_i+1) - sorted_R(idx_j);
        end
    end

    L_index = find(L_indicator > min(rate_common), 1) - 1;

    if isempty(L_index)
        % L_index = length(L_indicator);
        L_index = K - 1;
    end

    constants = zeros(K, 1);
    coeffMatrix = eye(K);

    if L_index == 0
        % disp("L_index is zero")
        C_portion = [min(rate_common) zeros(1, K-1)];
        C_portion = C_portion(reverse_idx);

    elseif L_index == K - 1

        for idx = 1:K-1
            coeffMatrix(idx,idx+1) = -1;
        end
        coeffMatrix(end,:) = 1;

        for idx = 1:K-1
            constants(idx) = sorted_R(idx + 1) - sorted_R(idx);
        end
        constants(end) = min(rate_common);

        solution = (coeffMatrix \ constants)';
        C_portion = solution(reverse_idx);

    else
        for idx = 1:L_index+1
            coeffMatrix(idx,idx+1) = -1;
        end
        coeffMatrix(end,:) = 1;

        for idx = 1:L_index+1
            constants(idx) = sorted_R(idx + 1) - sorted_R(idx);
        end

        constants(end) = min(rate_common);

        solution = (coeffMatrix \ constants)';
        C_portion = solution(reverse_idx);

    end

    C_portion(C_portion<0) = 0;

    % common message rate
    for idx_k = 1:para.K
        w_1(idx_k) = exp(-rate_common(idx_k)/alpha);
    end
    
    % common portion + private
    for idx_k = 1:para.K
        w_2(idx_k) = exp(-(rate_private(idx_k)+C_portion(idx_k))/alpha);
    end

    lagrangian_function.value{cnt_lag} = -alpha*log(sum(w_2))-lag_multi*(sum(C_portion)+alpha*log(sum(w_1)));
    % lagrangian_function.value{cnt_lag} = -alpha*log(sum(w_2));
    lagrangian_function.f{cnt_lag} = f_curr;
    lagrangian_function.C{cnt_lag} = C_portion;
    lagrangian_function.R{cnt_lag} = rate_private;
    lagrangian_function.common{cnt_lag} = rate_common;

end


[max_value, max_idx] = max(cell2mat(lagrangian_function.value));

% max_value;

f_GPI_RS = lagrangian_function.f{max_idx};
C_portion = lagrangian_function.C{max_idx};

end