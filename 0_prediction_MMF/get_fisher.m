function [fisher_matrix] = get_fisher(para, y_noisy, y_estimated, mu_est, nu_est, gain_est)

for i_ = 1:length(mu_est)
     l(i_).alpha = gain_est(i_);
     l(i_).mu = mu_est(i_);
     l(i_).tau = mu_est(i_) / para.scs;
     l(i_).nu = nu_est(i_);
end

fisher_matrix = zeros(length(mu_est)*4, length(mu_est)*4);
for i_ = 1:length(mu_est)
    for j_ = 1:length(mu_est)
        l1 = l(i_);
        l2 = l(j_);
        fisher_matrix(4*i_-3:4*i_, 4*j_-3:4*j_) = get_sub_matrix(para, l1, l2);
    end

end

function [sub_matrix] = get_sub_matrix(para, l1, l2)

sub_matrix = zeros(4, 4);
for m = 1:para.N
    for k = 1:para.Ns
        d1 = get_der(para, m, k, l1);
        d2 = get_der(para, m, k, l2);
        
        sub_matrix(1,1) = sub_matrix(1,1) + real(2*conj(d1.tau)*d2.tau+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.tau*d2.tau);
        sub_matrix(1,2) = sub_matrix(1,2) + real(2*conj(d1.tau)*d2.nu+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.tau*d2.nu);
        sub_matrix(1,3) = sub_matrix(1,3) + real(2*conj(d1.tau)*d2.g_real+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.tau*d2.g_real);
        sub_matrix(1,4) = sub_matrix(1,4) + real(2*conj(d1.tau)*d2.g_im+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.tau*d2.g_im);

        sub_matrix(2,1) = sub_matrix(2,1) + real(2*conj(d1.nu)*d2.tau+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.nu*d2.tau);
        sub_matrix(2,2) = sub_matrix(2,2) + real(2*conj(d1.nu)*d2.nu+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.nu*d2.nu);
        sub_matrix(2,3) = sub_matrix(2,3) + real(2*conj(d1.nu)*d2.g_real+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.nu*d2.g_real);
        sub_matrix(2,4) = sub_matrix(2,4) + real(2*conj(d1.nu)*d2.g_im+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.nu*d2.g_im);

        sub_matrix(3,1) = sub_matrix(3,1) + real(2*conj(d1.g_real)*d2.tau+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.g_real*d2.tau);
        sub_matrix(3,2) = sub_matrix(3,2) + real(2*conj(d1.g_real)*d2.nu+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.g_real*d2.nu);
        sub_matrix(3,3) = sub_matrix(3,3) + real(2*conj(d1.g_real)*d2.g_real+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.g_real*d2.g_real);
        sub_matrix(3,4) = sub_matrix(3,4) + real(2*conj(d1.g_real)*d2.g_im+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.g_real*d2.g_im);

        sub_matrix(4,1) = sub_matrix(4,1) + real(2*conj(d1.g_im)*d2.tau+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.g_im*d2.tau);
        sub_matrix(4,2) = sub_matrix(4,2) + real(2*conj(d1.g_im)*d2.nu+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.g_im*d2.nu);
        sub_matrix(4,3) = sub_matrix(4,3) + real(2*conj(d1.g_im)*d2.g_real+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.g_im*d2.g_real);
        sub_matrix(4,4) = sub_matrix(4,4) + real(2*conj(d1.g_im)*d2.g_im+conj(y_noisy(m, k)-y_estimated(m, k))*2*d1.g_im*d2.g_im);
    end
end


end

function d = get_der(para, m, k, l)

d.tau = l.alpha * exp(1j* 2*pi*(m-1)*l.nu)/sqrt(para.N) * (1j*2*pi*(k-1-(para.Ns-1)/2)*l.mu/l.tau) * exp(1j*2*pi*(k-1-(para.Ns-1)/2)*l.mu)*10^(para.SNR_per_element/20);
d.nu = l.alpha * 1j* 2*pi*(m-1) * exp(1j* 2*pi*(m-1)*l.nu)/sqrt(para.N) * exp(1j*2*pi*(k-1-(para.Ns-1)/2)*l.mu)*10^(para.SNR_per_element/20);
d.g_real = exp(1j* 2*pi*(m-1)*l.nu)/sqrt(para.N) * exp(1j*2*pi*(k-1-(para.Ns-1)/2)*l.mu)*10^(para.SNR_per_element/20);
d.g_im = 1j * exp(1j* 2*pi*(m-1)*l.nu)/sqrt(para.N) * exp(1j*2*pi*(k-1-(para.Ns-1)/2)*l.mu)*10^(para.SNR_per_element/20);

end


end