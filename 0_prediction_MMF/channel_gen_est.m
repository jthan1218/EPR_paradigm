function [H_ch,H_hat,H_hat_err_cov] = channel_gen_est(para)

H_ch = zeros(para.N,para.K);
H_hat = zeros(para.N,para.K);
H_hat_err_cov = zeros(para.N,para.N,para.K);

for idx_u = 1:para.K

    random_number_of_k = rand;

    %% true path information
    tmp = (randn(para.L,1) + 1i*randn(para.L,1))/sqrt(2);
    gain_true{idx_u} = para.sigma * (10.^(para.SNR/20)) .* tmp/norm(tmp);

    % gain_true{idx_u} = para.sigma * (10.^(10/10)) .* tmp/norm(tmp);
    %gain_true{idx_u} = gain.* exp(1j*unifrnd(0, 2*pi, para.L, 1));
    while true
        % nu_true{idx_u} = unifrnd(0.1, 0.15, [para.L, 1]);  %각도 좁힌거
        % nu_true{idx_u} = unifrnd(0, 1, [para.L, 1]);
        % nu_true{idx_u} = unifrnd(max(random_number_of_k-0.15,0), min(random_number_of_k+0.15, 0.75), [para.L, 1]);
        nu_true{idx_u} = unifrnd(max(random_number_of_k-0.05,0), min(random_number_of_k+0.05, 1.0), [para.L, 1]);
        % nu_true{idx_u}
        W_nu = nchoosek(nu_true{idx_u}, 2);
        D_nu = abs(W_nu(:,2)- W_nu(:, 1));
        min_inter_nu = min(D_nu);        
    
        while true
            
            mu_true{idx_u} = unifrnd(0, 1, [para.L, 1]); % defined as SCS times tau
            tau_true{idx_u} = mu_true{idx_u} / para.scs;
            W_mu = nchoosek(mu_true{idx_u}, 2);
            D_mu = abs(W_mu(:,2)- W_mu(:, 1));
            min_inter_mu = min(D_mu);
            if (min_inter_mu > 1/(2*para.Ns))
                %disp('possible mu found!!')
                break;
            elseif isempty(W_mu)
                break
            end
        
        end

        if (min_inter_nu > 1/para.N * 1/32) %각도 좁힌거
        % if (min_inter_nu) > 1/para.N * 1/2
            %disp('possible nu found!!')
            break;
        elseif isempty(W_nu)
            break;
        end
        
    end
    
    mu = mu_true{idx_u};
    nu = nu_true{idx_u};
    gain = gain_true{idx_u};

    path_sigma = sqrt(mean(norm(gain)^2));

    %% extract path information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mu, nu, gain, residue, y_noisy, y_estimated, y_full] = trivariate_NOMP(para, mu, nu, gain);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mu_est{idx_u} = mu;
    nu_est{idx_u} = nu;
    g = sqrt(path_sigma^2/2)*(randn(1)+1i*randn(1));
    gain_est_dl{idx_u} = para.eta * gain;
    gain_est_ul{idx_u} = gain;

    residue_est{idx_u} = residue;
    tau_est{idx_u} = (mu_est{idx_u}/para.scs);
    % same mu for DL and UL, lambda scaling for nu though.
    nu_est_dl{idx_u} = nu_est{idx_u}.*(para.ul_lda/para.dl_lda);
    nu_true_dl{idx_u} = nu_true{idx_u}.*(para.ul_lda/para.dl_lda);
    delta_f = para.dl_freq - para.ul_freq;
    

    %% reconstruct with estimated mu and nu & true mu and nu - pick a one subcarrier for representative
    % h_reconst_ul{idx_u} = zeros([para.N * para.Ns, 1]);
    h_reconst_dl{idx_u} = zeros([para.N * para.Ns, 1]);

    for l_ = 1: length(mu_est{idx_u})
        % h_reconst_ul{idx_u} = h_reconst_ul{idx_u} + kron(exp(1j*2*pi*(-(para.Ns-1)/2:(para.Ns-1)/2).' * mu_est{idx_u}(l_).')/sqrt(para.Ns), exp(1j*2*pi*(0:(para.N-1)).' * nu_est{idx_u}(l_).')/sqrt(para.N))* gain_est_ul{idx_u}(l_);
        h_reconst_dl{idx_u} = h_reconst_dl{idx_u} + exp(1i*2*pi*delta_f*tau_est{idx_u}(l_)).* kron(exp(1j*2*pi*(-(para.Ns-1)/2:(para.Ns-1)/2).' * mu_est{idx_u}(l_).')/sqrt(para.Ns), exp(1j*2*pi*(0:(para.N-1)).' * nu_est_dl{idx_u}(l_).')/sqrt(para.N))* gain_est_dl{idx_u}(l_);
    end

    % h_reconst_ul{idx_u} = reshape(h_reconst_ul{idx_u}, [para.N, para.Ns]);
    h_reconst_dl{idx_u} = reshape(h_reconst_dl{idx_u}, [para.N, para.Ns]);

    % pick a one subcarrier for representative
    % h_reconst_ul{idx_u} = h_reconst_ul{idx_u}(:, 1);
    h_reconst_dl{idx_u} = h_reconst_dl{idx_u}(:, 1);

    % h_true_ul{idx_u} = zeros([para.N * para.Ns, 1]);
    h_true_dl{idx_u} = zeros([para.N * para.Ns, 1]);

    for l_ = 1:para.L
        % h_true_ul{idx_u} = h_true_ul{idx_u} + kron(exp(1j*2*pi*(-(para.Ns-1)/2:(para.Ns-1)/2).' * mu_true{idx_u}(l_).')/sqrt(para.Ns), exp(1j*2*pi*(0:(para.N-1)).' * nu_true{idx_u}(l_).')/sqrt(para.N))* gain_true{idx_u}(l_);
        % h_true_dl{idx_u} = h_true_dl{idx_u} + exp(1i*2*pi*delta_f*tau_true{idx_u}(l_)).*  kron(exp(1j*2*pi*(-(para.Ns-1)/2:(para.Ns-1)/2).' * mu_true{idx_u}(l_).')/sqrt(para.Ns), exp(1j*2*pi*(0:(para.N-1)).' * nu_true_dl{idx_u}(l_).')/sqrt(para.N))* gain_true{idx_u}(l_);
        h_true_dl{idx_u} = h_true_dl{idx_u} + exp(1i*2*pi*delta_f*tau_true{idx_u}(l_)).*  kron(exp(1j*2*pi*(-(para.Ns-1)/2:(para.Ns-1)/2).' * mu_true{idx_u}(l_).')/sqrt(para.Ns), exp(1j*2*pi*(0:(para.N-1)).' * nu_true_dl{idx_u}(l_).')/sqrt(para.N))* (para.eta * gain_true{idx_u}(l_)+sqrt(1-para.eta^2)*g);
    end

    % h_true_ul{idx_u} = reshape(h_true_ul{idx_u}, [para.N, para.Ns]);
    h_true_dl{idx_u} = reshape(h_true_dl{idx_u}, [para.N, para.Ns]);

    % pick a one subcarrier for representative
    % h_true_ul{idx_u} = h_true_ul{idx_u}(:, 1);
    h_true_dl{idx_u} = h_true_dl{idx_u}(:, 1);

    %% H, est_H, Phi(w/ CRLB - fisher information, Jacobian)
    H_ch(:, idx_u) = h_true_dl{idx_u}/sqrt(10^(para.SNR_per_element/10));
    H_hat(:, idx_u) = h_reconst_dl{idx_u}/sqrt(10^(para.SNR_per_element/10));
    % H_hat_err_cov(:,:,idx_u) = diag(diag(get_G(para, mu_est{idx_u}, nu_est_dl{idx_u}, gain_est_ul{idx_u}, delta_f)'*pinv(get_fisher(para, y_noisy, y_estimated, mu_est{idx_u}, nu_est{idx_u}, gain_est_ul{idx_u}))*get_G(para, mu_est{idx_u}, nu_est_dl{idx_u}, gain_est_ul{idx_u}, delta_f)));
    % H_hat_err_cov(:,:,idx_u) = para.eta^2*diag(diag(get_G(para, mu_est{idx_u}, nu_est_dl{idx_u}, gain_est_dl{idx_u}, delta_f)'*pinv(get_fisher(para,y_noisy, y_estimated, mu_est{idx_u}, nu_est{idx_u}, gain_est_dl{idx_u}))*get_G(para, mu_est{idx_u}, nu_est_dl{idx_u}, gain_est_dl{idx_u}, delta_f)))+(1-para.eta^2)*eye(para.N);
    % H_hat_err_cov(:,:,idx_u) = para.eta^2*diag(diag(get_G(para, mu_est{idx_u}, nu_est_dl{idx_u}, gain_est_ul{idx_u}, delta_f)'*pinv(get_fisher(para,y_noisy, y_estimated, mu_est{idx_u}, nu_est{idx_u}, gain_est_ul{idx_u}))*get_G(para, mu_est{idx_u}, nu_est_dl{idx_u}, gain_est_ul{idx_u}, delta_f)))+(1-para.eta^2)*eye(para.N);
    H_hat_err_cov(:,:,idx_u) = para.eta^2*diag(diag(get_G(para, mu_est{idx_u}, nu_est{idx_u}, gain_est_ul{idx_u}, delta_f)'*pinv(get_fisher(para,y_noisy, y_estimated, mu_est{idx_u}, nu_est{idx_u}, gain_est_ul{idx_u}))*get_G(para, mu_est{idx_u}, nu_est{idx_u}, gain_est_ul{idx_u}, delta_f)))+(1-para.eta^2)*eye(para.N);
    H_hat_err_cov(:,:,idx_u) = H_hat_err_cov(:,:,idx_u)/para.N;
end



%% measure channel MSE characteristics
MSE = 0;
MSE_Phi = 0;
for idx_k = 1:para.K
    MSE = MSE + norm(H_ch(:, idx_k)-H_hat(:, idx_k))^2/norm(H_ch(:,idx_k))^2;
    MSE_Phi = MSE_Phi + trace(H_hat_err_cov(:,:,idx_k));
end
MSE = MSE/para.K;
MSE_Phi = MSE_Phi/para.K;

disp('------------w/o feedback-------------')
fprintf('generated channel MSE: %.2f\n', MSE);
fprintf('generated channel MSE w/ Phi: %.2f\n', MSE_Phi);
disp('-------------------------------------')
%%

end