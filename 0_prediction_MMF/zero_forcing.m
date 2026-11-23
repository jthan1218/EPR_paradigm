function f_ZF = zero_forcing(N,K,H_hat,SNR_N)

F = zeros(N,K);
for idx = 1:K
    F(:, idx) =  pinv(H_hat*H_hat'+ eye(N,N)*(1^2/SNR_N))*H_hat(:,idx); 
end

f_ZF = reshape(F,[N*K,1]);
f_ZF = f_ZF/norm(f_ZF,2);

end

