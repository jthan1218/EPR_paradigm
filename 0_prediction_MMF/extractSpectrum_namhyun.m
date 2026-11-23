function [mu_List, nu_List, gain_List, residue_List] = extractSpectrum_namhyun(y, para)

%% define 2D grid (mu, nu)
Ns = para.Ns;
Rs = para.OSR * para.Ns;
para.omega_mu = (0:(Rs-1))/Rs;

Na = para.N;
Ra = para.OSR * para.N;
para.omega_nu = (0:(Ra-1))/Ra;

mu_List = [];
nu_List = [];
gain_List = [];
y_r = y;

residue_List = [y_r' * y_r];

%% main loop - substract the most plausible argument in the grid
while true
    [max_arg_mu_new, max_arg_nu_new, max_gain_new, y_r, temp_res_info] = detectNew(y_r, para);
    if temp_res_info < para.tau
        break;
    end
    % fprintf('mu before refined: %.4f\n', max_arg_mu_new);
    % fprintf('nu before refined: %.4f\n', max_arg_nu_new);
    % fprintf('gain before refined: %.4f\n', max_gain_new);
    [max_arg_mu_new, max_arg_nu_new, max_gain_new, y_r] = refine0ne(y_r, max_arg_mu_new, max_arg_nu_new, max_gain_new);
    % fprintf('mu after refined: %.4f\n', max_arg_mu_new);
    % fprintf('nu after refined: %.4f\n', max_arg_nu_new);
    % fprintf('gain after refined: %.4f\n', max_gain_new);

    mu_List = [mu_List; max_arg_mu_new];
    nu_List = [nu_List; max_arg_nu_new];
    gain_List = [gain_List; max_gain_new];

    [mu_List, nu_List, gain_List, y_r] = refineAll(y_r, mu_List, nu_List, gain_List, para);


    residue_new = y_r'*y_r;
    residue_List = [residue_List; residue_new];
end

%% detect the most plausible one in the grid
function [max_arg_mu, max_arg_nu, max_gain, y_r, temp_res_info] = detectNew(y, para)
    i_ = 1;
    temp_power = 0;
    temp_res_info = 0;

    for mu = para.omega_mu
        j_ = 1;
        for nu = para.omega_nu
            mu_ = para.omega_mu(i_);
            nu_ = para.omega_nu(j_);
            u = kron(exp(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).' * mu_.')/sqrt(Ns), exp(1j* 2*pi*(0:(Na-1)).' * nu_.')/sqrt(Na));
            temp = abs(u'*y)^2;

            if temp > temp_power
                max_arg_mu = mu_;
                max_arg_nu = nu_;
                temp_power = temp;
                max_gain = u'*y;
            end
            
            if rem(i_, para.OSR) == 1 && rem(j_, para.OSR) == 1
                if temp > temp_res_info
                    temp_res_info = temp;
                end
                
            end

            j_ = j_ + 1;

        end
        i_ = i_ + 1;
    end

    y_r = y - max_gain * kron(exp(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).' * max_arg_mu')/sqrt(Ns), exp(1j* 2*pi*(0:(Na-1)).' * max_arg_nu')/sqrt(Na));
end

function [mu_new, nu_new, gain_new, y_r_next] = refine0ne(y_r, mu, nu, gain)
   dudm = kron(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).'.*exp(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).' * mu)/sqrt(Ns), exp(1j* 2*pi*(0:(Na-1)).' * nu)/sqrt(Na));
   d2udm = kron(-4*pi.^2*(-(Ns-1)/2:(Ns-1)/2).'.^2.*exp(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).' * mu)/sqrt(Ns), exp(1j* 2*pi*(0:(Na-1)).' * nu)/sqrt(Na));
   dudn = kron(exp(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).' * mu)/sqrt(Ns), 1j* 2*pi*(0:(Na-1)).'.*exp(1j* 2*pi*(0:(Na-1)).' * nu)/sqrt(Na));
   d2udn = kron(exp(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).' * mu)/sqrt(Ns), -4*pi.^2*(0:(Na-1)).'.^2.*exp(1j* 2*pi*(0:(Na-1)).' * nu)/sqrt(Na));
   d2udmdn = kron(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).'.*exp(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).' * mu)/sqrt(Ns), 1j* 2*pi*(0:(Na-1)).'.*exp(1j* 2*pi*(0:(Na-1)).' * nu)/sqrt(Na));
   
   der1 = [2*real(gain * y_r' * dudm); 2*real(gain * y_r' * dudn)];
   der2 = [2*real(gain * y_r' * d2udm) - 2 * abs(gain)^2*(dudm)'*(dudm), 2*real(gain * y_r' * d2udmdn) - 2 * abs(gain)^2*(dudn)'*(dudm); ... 
       2*real(gain * y_r' * d2udmdn) - 2 * abs(gain)^2*(dudm)'*(dudn), 2*real(gain * y_r' * d2udn) - 2 * abs(gain)^2*(dudn)'*(dudn)];
   y = y_r + gain * kron(exp(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).' * mu)/sqrt(Ns), exp(1j* 2*pi*(0:(Na-1)).' * nu)/sqrt(Na));
   der = pinv(der2)*der1;
   mu_new = real(mu - der(1));
   nu_new = real(nu - der(2));
   u_new = kron(exp(1j* 2*pi*(-(Ns-1)/2:(Ns-1)/2).' * mu_new')/sqrt(Ns), exp(1j* 2*pi*(0:(Na-1)).' * nu_new')/sqrt(Na));
   
   gain_new = u_new'*y;

   y_r_next = y - gain_new * u_new;

   if y_r_next'*y_r_next >= y_r'*y_r
    %fprintf('not effective! : before y_r: %.4f, but after y_r: %.4f \n', y_r_next'*y_r_next, y_r'*y_r);
    mu_new = mu;
    nu_new = nu;
    gain_new = gain;
   else
    %fprintf('refinement effective ! : before gain: %.4f but after gain: %.4f\n', gain, gain_new);
   end
    end

    function [mu_list, nu_list, gain_list, y_r] = refineAll(y_r, mu_list, nu_list, gain_list, para)
        K = length(mu_list);
        for i = 1:para.R_c
            order = 1:K;
            for j = 1:K
                l = order(j);
                mu = mu_list(l);
                nu = nu_list(l);
                gain = gain_list(l);
                [mu, nu, gain, y_r] = refine0ne(y_r, mu, nu, gain);
                mu_list(l) = mu;
                nu_list(l) = nu;
                gain_list(l) = gain;

            end
        end
    
    end

end