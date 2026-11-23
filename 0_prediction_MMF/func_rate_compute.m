function rate_matrix = func_rate_compute(ch, SNR_N, para, f_matrix)

K = para.K;
N = para.N;

rate_matrix = zeros(length(f_matrix), K);

for i = 1:length(f_matrix)

    if length(f_matrix{i}) == (K+1)*N
    
        % matrix computation
        for idx_k = 1:K
            h = ch(:, idx_k);
            phi = 0;

            A_k{idx_k}(:,:) = zeros(N*(K+1), N*(K+1));
            B_k{idx_k}(:,:) = zeros(N*(K+1), N*(K+1));
        
            for idx_k_inner = 1:1+K
                if idx_k_inner == 1
                    A_index = (idx_k_inner-1)*N + 1: (idx_k_inner-1)*N + N;
                    A_k{idx_k}(A_index, A_index) = zeros(N,N);
                else
                    A_index = (idx_k_inner-1)*N + 1: (idx_k_inner-1)*N + N;
                    A_k{idx_k}(A_index, A_index) = h*h' + phi;
                end
            end
        
            A_k{idx_k} = A_k{idx_k} + eye(N*(K+1), N*(K+1))*(1/SNR_N);
            B_k{idx_k} = A_k{idx_k};
            
            B_index = idx_k*N + 1: idx_k*N + N;
            B_k{idx_k}(B_index, B_index) = eye(N, N)*(1/SNR_N);
              
            C_k{idx_k}(:,:) = zeros(N*(K+1), N*(K+1));
            D_k{idx_k}(:,:) = zeros(N*(K+1), N*(K+1));
        
            for idx_k_inner = 1:1+K
                C_index = (idx_k_inner-1)*N + 1: (idx_k_inner-1)*N + N;
                C_k{idx_k}(C_index, C_index) = h*h';
            end
        
            C_k{idx_k} = C_k{idx_k} + eye(N*(K+1), N*(K+1))*(1^2/SNR_N);
            D_k{idx_k} = C_k{idx_k};
        
            D_index = 1:N;
            D_k{idx_k}(D_index, D_index) = eye(N, N)*(1^2/SNR_N);
        
        end
    
        for idx_k = 1:K
            SINR = real((f_matrix{i}'*A_k{idx_k}*f_matrix{i}/(f_matrix{i}'*B_k{idx_k}*f_matrix{i})))-1;
            rate_private(idx_k) = log2(1+SINR);
            rate_common(idx_k) = real(log2(f_matrix{i}'*C_k{idx_k}*f_matrix{i}/(f_matrix{i}'*D_k{idx_k}*f_matrix{i})));
        end
    
        R = rate_private;
    
    else % no RS
        
        for idx_k = 1:K
            h = ch(:, idx_k);
            phi = 0;
            
            A{idx_k}(:,:) = zeros(N*K, N*K);
            B{idx_k}(:,:) = zeros(N*K, N*K);
        
            for idx_k_inner = 1:K
                A_index = (idx_k_inner-1)*N + 1: (idx_k_inner-1)*N + N;
                A{idx_k}(A_index, A_index) = h*h' + phi;
            end
        
            A{idx_k} = A{idx_k} + eye(N*K, N*K)*(1^2/SNR_N);
            B{idx_k} = A{idx_k};
            
            B_index = (idx_k-1)*N + 1: (idx_k-1)*N + N;
            B{idx_k}(B_index, B_index) = eye(N, N)*(1^2/SNR_N);

        end
    
        for idx_k = 1:K
            SINR = real(f_matrix{i}'*A{idx_k}*f_matrix{i}/(f_matrix{i}'*B{idx_k}*f_matrix{i}))-1;
            rate(idx_k) = log2(1+SINR);
        end
    
        R = rate;
        
    end

    rate_matrix(i, :) = R;

end
    
end