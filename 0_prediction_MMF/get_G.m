function G = get_G(para, mu_est, nu_est, gain_est, f)

for i_ = 1:length(mu_est)
     l(i_).alpha = gain_est(i_);
     l(i_).mu = mu_est(i_);
     l(i_).tau = mu_est(i_) / para.scs;
     l(i_).nu = nu_est(i_);
end

G = zeros(length(mu_est)*4, para.N);
for i_ = 1:para.N
    for j_ = 1:length(mu_est)
        G(4*j_-3:4*j_, i_) = [1j*2*pi*f * l(j_).alpha * exp(1j*2*pi*f*l(j_).tau) * exp(1j*2*pi*l(j_).nu*(i_-1))/sqrt(para.N); 
            l(j_).alpha * 1j*2*pi*(i_-1) * exp(1j*2*pi*(i_-1)*l(j_).nu)/sqrt(para.N) * exp(1j*2*pi*f*l(j_).tau);
            exp(1j*2*pi*f*l(j_).tau) * exp(1j*2*pi*l(j_).nu*(i_-1))/sqrt(para.N);
            1j * exp(1j*2*pi*f*l(j_).tau) * exp(1j*2*pi*l(j_).nu*(i_-1))/sqrt(para.N)];
    end
    
end