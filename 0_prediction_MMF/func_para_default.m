function para = func_para_default

para.N = 6;
para.K = 3;

para.L = 5;

para.iter = 30;
para.P_db_Q = 10:5:40;
para.sigma = 1;

para.lagr_max = 5;
para.lagr_step = 0.25;
para.th = 0.01;
para.mmse_th = 0.1;
para.alpha = 0.1;
para.Branch = 'Rate_Splitting:Branch_1';

para.SAF_M = 10;
para.mmse_th = 0.1;

para.alpha_increase = 0.1;
para.loop_terminate = 50;


para.Ns = 128; % number of subcarriers
para.OSR = 2; % over sampling rate
para.R_c = 3;
% para.scs = 400 * 1.0e+03; % 30kHz, https://www.sharetechnote.com/html/5G/5G_FrameStructure.html#RadioFrame_Structure_SlotConfig0_u_0
para.scs = 400 * 1.0e+03;
para.symbol_duration = 1/para.scs;
para.cp_duration = para.symbol_duration;%2.34 * 1.0e-06;
para.c = 299792358; % in [m/s] original
% para.c = 300000000;
% para.ul_freq = 7.125e+09; %original
% para.ul_freq = 7.15e+09;
para.ul_freq = 2.535e+09;
para.ul_lda = para.c/para.ul_freq;
para.d = (para.c/para.ul_freq)/2;
% para.dl_freq = para.ul_freq+0.1e+9; % original
para.dl_freq = para.ul_freq+0.12e+09;
% para.dl_freq = para.ul_freq+0.15e+9;
para.dl_lda = para.c/para.dl_freq;
para.SNR_per_element = 10;
para.SNR = para.SNR_per_element + 10*log10(para.N*para.Ns);
dl_freqs = linspace(para.ul_freq, para.ul_freq+10e+08, 51);
dl_freqs = dl_freqs(2: length(dl_freqs));

end