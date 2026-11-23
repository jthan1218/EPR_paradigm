function f_noRS = without_RS(N,K,H_hat,H_hat_err_cov,SNR_N, para)

%% without RS

% %% input ZF version
% F_fd_MSE = pinv(H_hat*H_hat' + (1/SNR_N)*eye(N, N))*H_hat;
% 
% f_curr = reshape(F_fd_MSE,N*K,1);
% f_curr = f_curr/norm(f_curr, 2);

%% input MRT version
f_mtx = zeros(N,K);

for idx_k = 1:K
    f_mtx(:, idx_k) = H_hat(:, idx_k);
end

f = reshape(f_mtx, [N*K, 1]);
f_curr = f/norm(f, 2);

%%
rate_noRS = zeros(1, K);

for idx_k = 1:K
    h = H_hat(:, idx_k);
    phi = H_hat_err_cov(:, :, idx_k);
    % phi = 0;
    
    A{idx_k}(:,:) = zeros(N*K, N*K);
    B{idx_k}(:,:) = zeros(N*K, N*K);

    for idx_k_inner = 1:K
        A_index = (idx_k_inner-1)*N + 1: (idx_k_inner-1)*N + N;
        A{idx_k}(A_index, A_index) = h*h' + phi;
    end

    A{idx_k} = A{idx_k} + eye(N*K, N*K)*(1^2/SNR_N);
    B{idx_k} = A{idx_k};
    
    B_index = (idx_k-1)*N + 1: (idx_k-1)*N + N;
    B{idx_k}(B_index, B_index) = eye(N, N)*(1^2/SNR_N) + phi;    

end

t = 0;
alpha = para.alpha;

while 1
    t = t + 1;
    f_prev = f_curr;
    f = f_curr;

    for idx_k = 1:K
        SINR = real(f'*A{idx_k}*f/(f'*B{idx_k}*f)) - 1;
        rate_noRS(idx_k) = log2(1+SINR);
        w(idx_k) = exp(-(rate_noRS(idx_k))/alpha);
    end

    w = w/sum(w);

    X_kkt = 0;
    Y_kkt = 0;

    for idx_k = 1:para.K
        X_inner{idx_k} = w(idx_k)*A{idx_k}/(f'*A{idx_k}*f);
        Y_inner{idx_k} = w(idx_k)*B{idx_k}/(f'*B{idx_k}*f);
        X_kkt = X_kkt + X_inner{idx_k};
        Y_kkt = Y_kkt + Y_inner{idx_k};
    end

    for idx_k = 1:K
        Y_kkt_part{idx_k} = Y_kkt((idx_k-1)*N + 1:(idx_k)*N,(idx_k-1)*N + 1:(idx_k)*N);
        inv_Y_kkt_part{idx_k} = pinv(Y_kkt_part{idx_k});
    end

    inv_Y_kkt = blkdiag(inv_Y_kkt_part{:});

    f_curr = pinv(Y_kkt)*X_kkt*f;
    % f_curr = inv_Y_kkt * X_kkt * f;
    f_curr = f_curr/norm(f_curr, 2);

    if norm(f_prev-f_curr, 2) < para.th
        break;
    end

    if t > para.loop_terminate * 2
        alpha = alpha*(1+para.alpha_increase);
        % para.th = para.th * 1.01;
        % para.th = para.th * 1.05;
        para.th = para.th * 1.01; % (-0.5, 0.5)에서는 없이, (0, 1)에서는 설정 %rate 뽑을때 설정
        % para.th = para.th * 1.75; % latency 뽑을 때 설정

        t = 0;
        f = randn(para.N*para.K,1) + 1i*randn(para.N*para.K,1);

        f = f/norm(f,2);

        f_curr = f;
    end

end

f_noRS = f_curr;

end