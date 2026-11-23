function f_MRT = MRT(N,K, H_hat)

F = zeros(N,K);

for idx_k = 1:K
    F(:, idx_k) = H_hat(:, idx_k);
end

f_MRT = reshape(F, [N*K, 1]);
f_MRT = f_MRT/norm(f_MRT, 2);

end