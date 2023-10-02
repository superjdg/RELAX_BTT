%20221124
%@article{wang2023min, 
% title={Min-max Probe Placement and Extended Relaxation Estimation Method for Processing Blade Tip Timing Signals}, 
% author={Wang, Zengkun and Petre, Stoica and Dave, Zachariah and Prabhu, Babu and Zhibo, Yang}, 
% journal={IEEE TRANSACTIONS ON INSTRUMENTATION AND MEASUREMENT}, 
% year={2023}}
clear all
close all

%parameter setting
rng(666)
omega = 6000/60; %rotating speed
V = 24;%virtual probe number
P = 4;
Q = 60;                %number of revs
M = Q*P;       %before padding (actual samples)
N = Q*V;
delta_t = 1/omega/V; %virtual sampling interval

%signal generation
K = 4;
probe_lay = [0 1 3 10];   %probe layout
f = sort(randperm(V/2-1,K)'*omega+randperm(omega/4,K)'*3+normrnd(0,1,K,1));%3 is not common divisor of rotating fre
f(3) = f(2)+1;

A = rand(length(f),1)*0.2+1;                  %a_k
A(1) = rand(1,1)*0.2+0.3;
A(4) = rand(1,1)*0.2+0.4;
phase = rand(length(f),1)*pi/2-pi/4;              %phase
t = zeros(Q*P,1);
for i = 0:Q-1
    for j = 1:length(probe_lay)
        t(i*P+j) = (V*i+probe_lay(j))*delta_t; %the index of the actual sampling times
    end
end
x_pure = zeros(size(t));
for i = 1:length(f)
    x_pure=x_pure+A(i)*sin(2*pi*f(i)*t+phase(i)); %virtual displacement
end
power_sig = mean(x_pure.^2);

SNRs = [0:5:25];
N_mc =504;
MSEs = zeros(length(SNRs)*K,3*4);
Results = zeros(length(SNRs)*N_mc,K*4);
for i_snr = 1:length(SNRs)
    SNR = SNRs(i_snr);
    waitbar(i_snr/length(SNRs))
    var_noi = power_sig/10^(SNR/10);
    MSE_relax = zeros(K,3);
    MSE_BOMP = MSE_relax;
    MSE_music = MSE_relax;
    tic
    parfor n_mc = 1:N_mc
%         disp(n_mc)
        seed = n_mc+(N_mc-1)*SNR;
        rng(2022+seed)

        x = zeros(length(x_pure),1);
        noise = normrnd(0,sqrt(var_noi),length(x),1);
        x = x_pure+noise;
%         tic
        %RELAX
        [Thitas,Var_s,f_fft] = RELAX_BTT_no_syn(x,t,probe_lay,V,K,round(1000000/V)*V);
        Thitas_cal = zeros(3*K,K);%FRE, AMP, PHASE
        for k = 1:K
            Thitas_cal(3*k-2,1:k) = Thitas(3*k-2,1:k);
            Thitas_cal(3*k-1,1:k) = sqrt(Thitas(3*k-1,1:k).^2+Thitas(3*k,1:k).^2);
            Thitas_cal(3*k,1:k) = atan((Thitas(3*k,1:k)./Thitas(3*k-1,1:k)));
        end
        [f_relax_sort,temp_index] = sort(Thitas_cal(3*K-2,1:K));
        A_relax_sort = Thitas_cal(3*K-1,temp_index);
        Phase_relax_sort = Thitas_cal(3*K,temp_index);
        f_err = (f_relax_sort-f(1:K)').^2;
        A_err = (A_relax_sort-A(1:K)').^2;
        Phase_err = (Phase_relax_sort-phase(1:K)').^2;
        MSE_relax = MSE_relax+[f_err' A_err' Phase_err']./N_mc;
%         Results((i_snr-1)*N_mc+n_mc,1:K) = [f_relax_sort];
%        toc
        %Block-OMP
        OMP_A = zeros(length(t),length(f_fft)); %Ö»ËãÒ»°ëµÄÆµÆ×
        for i = 1:length(f_fft)/2
            OMP_A(:,2*i-1) = sin(2*pi*f_fft(i)*t);%sin
            OMP_A(:,2*i) = cos(2*pi*f_fft(i)*t);%cos
        end
        [thita_OMP,residual_OMP] = Block_OMP( OMP_A, x,length(f),10^(-3));
        cor_index = find(thita_OMP);
%         figure()
%         plot(f_fft(1:length(f_fft)/2),sqrt(thita_OMP(1:2:end).^2+thita_OMP(2:2:end).^2))
        Phase_omp = zeros(1,K);
        A_omp = Phase_omp;
        f_omp = Phase_omp;
        Phase_omp_temp = zeros(1,length(f));
        A_omp_temp = Phase_omp_temp;
        f_omp_temp = Phase_omp_temp;
        for i = 1:length(cor_index)/2
            f_omp_temp(i) = f_fft(cor_index(2*i)/2);
            A_omp_temp(i) = sqrt(sum(thita_OMP(cor_index(2*i-1:2*i)).^2));
            Phase_omp_temp(i) =atan((thita_OMP(cor_index(2*i))./thita_OMP(cor_index(2*i-1))));
        end
        [f_omp_sort,temp_index] = sort(f_omp_temp);
        A_omp_sort = A_omp_temp(temp_index);
        Phase_omp_sort = Phase_omp_temp(temp_index);
        f_err = (f_omp_sort-f(1:K)').^2;
        A_err = (A_omp_sort-A(1:K)').^2;
        Phase_err = (Phase_omp_sort-phase(1:K)').^2;
        MSE_BOMP = MSE_BOMP+[f_err' A_err' Phase_err']./N_mc;
%         Results((i_snr-1)*N_mc+n_mc,K+1:2*K) = [f_omp_sort];
%         toc
        %MUSIC
        m = Q/2*P;
        S = zeros(m,Q/2+1);
        for i = 0:Q/2
            S(:,i+1) = x(i*P+1:i*P+m);%[x(M-i*P+1:M);x(1:M-i*P)];
        end
        [spectrum] = MUSIC_sparse(S,t(1:m),f_fft(1:length(f_fft)/2), K);
        [pks,idx1] = findpeaks(-spectrum);
        figure()
        plot(f_fft(1:length(f_fft)/2),spectrum)
        f_pks = f_fft(idx1);
        [maxs,idx2] = maxk(pks,K);
        f_music = f_pks(idx2);
        [f_music_sort,temp_index] = sort(f_music);
        f_err = (f_music_sort-f(1:K)').^2;
        A_err = zeros(size(f_err));
        Phase_err = A_err;
        MSE_music = MSE_music+[f_err' A_err' Phase_err']./N_mc;
%         Results((i_snr-1)*N_mc+n_mc,2*K+1:3*K) = [f_music_sort];
%         Results((i_snr-1)*N_mc+n_mc,3*K+1:4*K) = [f(1:K)'];
%         toc
%         toc
    end
    toc
    %CRB
    derivative_matrix = zeros(M,3*length(f));
    for i = 1:length(f)
        %for frequency
        derivative_matrix(:,(i-1)*3+1) = 2*pi*t*A(i).*cos(2*pi*f(i).*t+phase(i));
        %Derivation of amplitude
        derivative_matrix(:,(i-1)*3+2) = sin(2*pi*f(i)*t+phase(i));
        %for phase
        derivative_matrix(:,(i-1)*3+3) = A(i)*cos(2*pi*f(i)*t+phase(i));
    end
    FIM = derivative_matrix'*derivative_matrix/var_noi;
    CRB_vector = diag(inv(FIM));
    CRB_f =  CRB_vector(1:3:3*K-2);
    CRB_a =  CRB_vector(2:3:3*K-1);
    CRB_phase=  CRB_vector(3:3:3*K);
    var_crb = [CRB_f CRB_a CRB_phase];
    MSEs((i_snr-1)*K+1:(i_snr-1)*K+K,1:3) = MSE_relax;
    MSEs((i_snr-1)*K+1:(i_snr-1)*K+K,4:6) = MSE_BOMP;
    MSEs((i_snr-1)*K+1:(i_snr-1)*K+K,7:9) = MSE_music;
    MSEs((i_snr-1)*K+1:(i_snr-1)*K+K,10:12) = var_crb;
end
save('MSEs_snr_white.mat','f','A','phase','SNRs','MSEs')





    
    