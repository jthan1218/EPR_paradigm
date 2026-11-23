function [mu_est, nu_est, gain_est, residue, y_noisy, y_estimated, y_full] = trivariate_NOMP(para, mu_true, nu_true, gain_true)

set(0,'DefaultLineMarkerSize',12);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',14);

% original signal
y_full = zeros(para.Ns*para.N, 1);
for k_ = 1:para.L
    y_full = y_full + kron(exp(1j* 2*pi*(-(para.Ns-1)/2:(para.Ns-1)/2).' * mu_true(k_).')/sqrt(para.Ns), exp(1j* 2*pi*(0:(para.N-1)).' * nu_true(k_).')/sqrt(para.N)) * gain_true(k_);
end
%% Choose measurement type
% normal measurements
M = para.Ns*para.N; % number of measurements = N
S = eye(M);

%% Measurements
% add measurement noise
noise = para.sigma*(randn(M,1) + 1j*randn(M,1))/sqrt(2);
y_noisy = y_full + noise;

% take measurements
y = S * y_noisy;
    
%% Algorithm 
% stopping criterion parameter tau
% when we know the noise level sigma^2,
% we choose tau to ensure that the false alarm rate is fixed:
% ie., the probability that we over estimate the size of the
% support set is fixed to some pre-determined number p_fa
para.p_fa = 1e-2;
para.tau = para.sigma^2 * ( log(M) - log( log(1/(1-para.p_fa)) ) );

% fprintf('max mu : %.4f\n', mu_true(1));
% fprintf('max nu : %.4f\n', nu_true(1));
% fprintf('max gain : %.4f\n', gain_true(1));

[mu_est, nu_est, gain_est, residue] = extractSpectrum_namhyun(y, para);

% original signal
y_estimated = zeros(para.Ns*para.N, 1);
for k_ = 1:length(mu_est)
    y_estimated = y_estimated + kron(exp(1j* 2*pi*(-(para.Ns-1)/2:(para.Ns-1)/2).' * mu_est(k_).')/sqrt(para.Ns), exp(1j* 2*pi*(0:(para.N-1)).' * nu_est(k_).')/sqrt(para.N)) * gain_est(k_);
end

y_estimated = reshape(y_estimated, [para.N, para.Ns]);
y_noisy = reshape(y_estimated, [para.N, para.Ns]);
y_full = reshape(y_full, [para.N, para.Ns]);
% %% Plot results
% figure; 
% subplot(1,3,1); polar(mu_true, abs(gain_true), 'bo');
% hold on; polar(mu_est, abs(gain_est),'rx');
% title(sprintf('Magnitude and frequency\nof each sinusoid'));
% legend({'Truth', 'Estimate'}, 'Location', 'SouthOutside');
% 
% subplot(1,3,2); polar(nu_true, abs(gain_true), 'bo');
% hold on; polar(nu_est, abs(gain_est),'rx');
% title(sprintf('Magnitude and frequency\nof each sinusoid'));
% legend({'Truth', 'Estimate'}, 'Location', 'SouthOutside');
% 
% subplot(1,3,3); polar(angle(gain_true), abs(gain_true), 'bo');
% hold on; polar(angle(gain_est), abs(gain_est), 'rx');
% title(sprintf('Magnitude and phase\nof each sinusoid'));
% legend({'Truth', 'Estimate'}, 'Location', 'SouthOutside');
% 
% 
% figure; stem(0:(length(residue)-1), 10*log10(residue/para.sigma^2/M));
% xlabel('Number of detected sinusoids');
% ylabel('Residual power (relative to \sigma^2) in dB');
% hold on; plot(para.L, 0, 'rx');
% legend({'Estimate', 'Truth'});
