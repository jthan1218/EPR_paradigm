clear all;
close all;
warning 'off'
clc;

%% ==================== [자동 리소스 프로파일링] 시작 ====================
% 목적: 이 스크립트의 (1) 총 실행 시간, (2) Peak RAM(프로세스 RSS) 자동 측정
% - macOS/Linux: `ps -o rss=` , Windows: `tasklist` 로 프로세스 실제 점유 메모리 샘플링
% - 백그라운드 timer 가 0.5초마다 RSS 를 읽어 최댓값을 setappdata(root)에 누적
% - clear all 이후에 배치해야 타이머 변수가 지워지지 않음 (appdata 는 clear all 에 영향 없음)
delete(timerfindall('Tag','prof_mem'));          % 이전 실행에서 남은 타이머 정리
prof_t_start = tic;                               % 실행 시간 측정 시작
prof_pid     = feature('getpid');                % 현재 MATLAB 프로세스 PID
setappdata(0, 'prof_peak_kb', 0);                % Peak RSS(KB) 저장소 초기화
prof_timer = timer('Tag','prof_mem', ...
    'ExecutionMode','fixedRate', 'Period',0.5, ...
    'TimerFcn', @(~,~) prof_sample_mem(prof_pid));
start(prof_timer);                               % 메모리 샘플링 백그라운드 시작
%% ======================================================================

%% Parameters

%% true channel인지, channel reconstruction인지 구별하고, channel reconstruction이면 ECM 까먹지 말것.
para = func_para_default;
para.P_db_Q = 0:5:20;
para.iter = 15;

para.N = 12;
para.K = 4;

N = para.N;
K = para.K;

para.eta = 0.95;

%% Minimum rate parameters
iter = para.iter;
% lag_multi_range = 0.2:0.05:0.4;
% lag_multi_range = 0.1:0.02:0.4;  % 현재 reference, 근데 원래 MMF 논문에서도 minimum value는 0이네
% lag_multi_range = 0.1:0.08:0.5; % Not bad
% lag_multi_range = 0.15:0.05:0.6; % 여기서 high SNR 갈수록 lag value 낮아진다는 것 확인
% lag_multi_range = 0.4:0.05:1.0;
lag_multi_range = 0.5:0.1:0.9;
% lag_multi_range = 0.1:0.05:0.4;  % 안테나 많아지면 너무 오래걸려서 step 증가시킴.(N, K = 32, 8 Case)
% lag_multi_range = 0.2:0.02:0.36; %WMMSE랑 거의 비슷, proposed method가 살짝 outperform
% lag_multi_range = 0.01:0.01:0.4; %WMMSE랑 거의 비슷, proposed method가 살짝 outperform
SNR_dB_range = para.P_db_Q;

% 왜 lag_multi_range를 1까지 설정하면 값이 이상해지나, lagrangian_function으로 비교해서 큰 값만 뽑아오면 되는데...

%% Figure parameters
min_rate_MRT = zeros(iter,length(SNR_dB_range));
min_rate_zf = zeros(iter,length(SNR_dB_range));
min_rate_wo = zeros(iter,length(SNR_dB_range));
MMF_min_sumrate_RS = zeros(iter,length(SNR_dB_range));
min_sumrate_WMMSE = zeros(iter,length(SNR_dB_range));
MMF_min_sumrate_RS_noECM = zeros(iter,length(SNR_dB_range));
Time_GPI = zeros(iter,length(SNR_dB_range));
Time_WMMSE = zeros(iter,length(SNR_dB_range));
Time_noRS = zeros(iter, length(SNR_dB_range));

%% The minimum rate
cnt_SNR = 0;
for SNR_dB = SNR_dB_range
    cnt_SNR = cnt_SNR + 1;
    SNR_N = 10^(SNR_dB/10);

    status = cnt_SNR/length(SNR_dB_range);

    disp(status);

    for iter=1:para.iter

        %% Channel generation
        [H_ch,H_hat,H_hat_err_cov] = channel_gen_est(para);

        %% MRT
        f_MRT = MRT(N,K,H_hat);

        %% Zero forcing
        f_ZF = zero_forcing(N,K,H_hat, SNR_N);

        %% WMMSE
        Time_start.wmmse = cputime;
        [f_WMMSE, C_WMMSE] = WMMSE_ver6(SNR_dB,K,N,H_hat,para);
        Time_end.wmmse = cputime;
        Time.wmmse = Time_end.wmmse - Time_start.wmmse;

        %% Without-RS_ch_est
        Time_start.noRS = cputime;
        f_noRS = without_RS(N,K,H_hat,H_hat_err_cov,SNR_N, para);
        Time_end.noRS = cputime;
        Time.noRS = Time_end.noRS - Time_start.noRS;

        %% MMF_GPI_RSMA
        Time_start.proposed = cputime;
        [f_GPI_RS, C_GPI_RSMA, max_idx] = MMF_GPI_RS(lag_multi_range,K,N,H_hat,H_hat_err_cov,SNR_N, para); %perfect CSI로 변경
        lag_value = lag_multi_range(max_idx)
        Time_end.proposed = cputime;
        Time.proposed = Time_end.proposed - Time_start.proposed;

        %% MMF_GPI_RSMA_noECM
        Time_start.noECM = cputime;
        [f_GPI_RS_noECM, C_GPI_RSMA_noECM, max_idx_2] = MMF_GPI_RS_noECM(lag_multi_range,K,N,H_hat,H_hat_err_cov,SNR_N, para); %perfect CSI로 변경
        lag_value_noECM = lag_multi_range(max_idx_2)
        Time_end.noECM = cputime;

        Time.noECM = Time_end.noECM - Time_start.noECM;
        %% Rate computation
        f_matrix{1} = f_MRT;
        f_matrix{2} = f_ZF;
        f_matrix{3} = f_noRS;
        f_matrix{4} = f_GPI_RS;
        f_matrix{5} = f_WMMSE;
        f_matrix{6} = f_GPI_RS_noECM;

        rate_matrix = func_rate_compute(H_ch, SNR_N, para, f_matrix);

        min_rate_MRT(iter, cnt_SNR) = min(rate_matrix(1,:));
        min_rate_zf(iter, cnt_SNR) = min(rate_matrix(2,:));
        min_rate_wo(iter, cnt_SNR) = min(rate_matrix(3,:));
        rate_check = rate_matrix(3,:);
        % rate_check = rate_matrix(4,:)+C_GPI_RSMA
        MMF_min_sumrate_RS(iter, cnt_SNR) = min(rate_matrix(4,:)+C_GPI_RSMA)
        min_sumrate_WMMSE(iter, cnt_SNR) = min(rate_matrix(5,:)+C_WMMSE')
        MMF_min_sumrate_RS_noECM(iter, cnt_SNR) = min(rate_matrix(6,:)+C_GPI_RSMA_noECM)

        MMF_check = mean(MMF_min_sumrate_RS(:,cnt_SNR))
        MMF_check_WMMSE = mean(min_sumrate_WMMSE(:,cnt_SNR))
        MMF_check_noECM = mean(MMF_min_sumrate_RS_noECM(:,cnt_SNR))

        Time_GPI(iter, cnt_SNR) = Time.proposed;
        Time_WMMSE(iter, cnt_SNR) = Time.wmmse;
        Time_noRS(iter, cnt_SNR) = Time.noRS;

    end
end

mean_min_rate_MRT = mean(min_rate_MRT);
mean_min_rate_zf = mean(min_rate_zf);
mean_min_rate_WMMSE = mean(min_sumrate_WMMSE);
mean_min_rate_wo= mean(min_rate_wo);
mean_MMF_RS = mean(MMF_min_sumrate_RS);
mean_MMF_RS_noECM = mean(MMF_min_sumrate_RS_noECM);

Time_avg.proposed = mean(Time_GPI);
Time_avg.wmmse = mean(Time_WMMSE);
Time_avg.noRS = mean(Time_noRS);

computation_time_proposed = mean(Time_avg.proposed);
computation_time_wmmse = mean(Time_avg.wmmse);
computation_time_noRS = mean(Time_avg.noRS);

%% Figure
figure(1)
p1 = plot(SNR_dB_range,mean_MMF_RS,'-o','color',[1.00,0.00,0.00],'LineWidth',3, 'MarkerSize', 14);
hold on;
grid on;
p2 = plot(SNR_dB_range,mean_min_rate_WMMSE,'-*','color',[0.5,0.5,1.00],'LineWidth',3, 'MarkerSize', 14);
p3 = plot(SNR_dB_range,mean_MMF_RS_noECM,'--o','color',[1,0,0],'LineWidth',3, 'MarkerSize', 14);
p4 = plot(SNR_dB_range,mean_min_rate_wo,'-x','color',[0,0.75,1.0],'LineWidth',3, 'MarkerSize', 14);
p5 = plot(SNR_dB_range,mean_min_rate_zf,'-d','color',[0.00,0.00,0.00],'LineWidth',3, 'MarkerSize', 14);
p6 = plot(SNR_dB_range,mean_min_rate_MRT,'-<','color',[0.50,0.00,0.50],'LineWidth',3, 'MarkerSize', 14);

% 플롯 순서 재정렬: MRT, RZF, Without RS, Proposed (w/o ECM), WMMSE, Proposed
uistack(p6, 'top'); % MRT
uistack(p5, 'top'); % RZF
uistack(p4, 'top'); % Without RS
uistack(p3, 'top'); % Proposed (w/o ECM)
uistack(p2, 'top'); % WMMSE
uistack(p1, 'top'); % Proposed

% Legend를 원래 순서대로
legend([p1, p2, p3, p4, p5, p6], 'Proposed', 'WMMSE', 'Proposed (w/o ECM)', 'GPI (w/o RS)', 'RZF', 'MRT', 'Location', 'northwest');



% legend('Proposed','Without RS','RZF','MRT','WMMSE','Location','northwest');
% legend('Proposed','WMMSE','Proposed (w/o ECM)','GPI (w/o RS)','RZF','MRT','Location','northwest');
xlabel('SNR (dB)');
ylabel('The Minimum Spectral Efficiency (bps/Hz)');

% figure(2)
% plot(SNR_dB_range, Time_avg.proposed, '-o','color',[1.00,0.00,0.00],'LineWidth',1.2);
% hold on;
% plot(SNR_dB_range, Time_avg.wmmse,'-*','color',[0.50,0.50,1.00],'LineWidth',1.2);
% legend('Proposed MMF-GPI-RS','WMMSE');
set(gca, 'FontSize', 16);

%% ==================== [자동 리소스 프로파일링] 종료 ====================
prof_elapsed = toc(prof_t_start);                % 총 실행 시간 (초)
stop(prof_timer); delete(prof_timer);            % 백그라운드 타이머 정리
prof_peak_kb = getappdata(0, 'prof_peak_kb');    % 측정된 Peak RSS (KB)
prof_peak_MB = prof_peak_kb / 1024;
prof_peak_GB = prof_peak_MB / 1024;

fprintf('\n');
fprintf('========================================================\n');
fprintf('            [실행 리소스 프로파일링 결과]\n');
fprintf('========================================================\n');
fprintf('  총 실행 시간 (Elapsed Time) : %.2f 초  (%.2f 분)\n', prof_elapsed, prof_elapsed/60);
fprintf('  Peak RAM  (Peak RSS)       : %.1f MB  (%.3f GB)\n', prof_peak_MB, prof_peak_GB);
fprintf('  측정 프로세스 PID           : %d\n', prof_pid);
fprintf('========================================================\n\n');

% result.mat 오염 방지를 위해 프로파일링 임시 변수 정리
clear prof_timer prof_pid prof_t_start prof_peak_kb;
%% ======================================================================

save('result')


%% ---- [프로파일링 보조 함수] 프로세스 RSS(KB) 샘플링 ----
% 본 함수는 스크립트 로컬 함수이므로 파일 맨 끝에 위치해야 함 (R2016b 이상)
function prof_sample_mem(pid)
    kb = NaN;
    try
        if ispc
            % Windows: tasklist 출력에서 메모리 사용량(K) 추출
            [st, out] = system(sprintf('tasklist /FI "PID eq %d" /FO CSV /NH', pid));
            if st == 0
                tok = regexp(out, '"([\d.,]+) K"', 'tokens');
                if ~isempty(tok)
                    kb = str2double(regexprep(tok{end}{1}, '[.,]', ''));
                end
            end
        else
            % macOS / Linux: ps 의 RSS(KB) 직접 사용
            [st, out] = system(sprintf('ps -o rss= -p %d', pid));
            if st == 0
                kb = str2double(strtrim(out));
            end
        end
    catch
        return;   % 샘플링 실패는 무시 (본 계산에 영향 주지 않도록)
    end

    if ~isnan(kb) && kb > getappdata(0, 'prof_peak_kb')
        setappdata(0, 'prof_peak_kb', kb);
    end
end
