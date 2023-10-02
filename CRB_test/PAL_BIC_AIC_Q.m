%20220110

clear all
close all
%@article{wang2023min, 
% title={Min-max Probe Placement and Extended Relaxation Estimation Method for Processing Blade Tip Timing Signals}, 
% author={Wang, Zengkun and Petre, Stoica and Dave, Zachariah and Prabhu, Babu and Zhibo, Yang}, 
% journal={IEEE TRANSACTIONS ON INSTRUMENTATION AND MEASUREMENT}, 
% year={2023}}
%parameter setting
rng(666)
omega = 6000/60; %rotating speed
V = 24;%virtual probe number
probe_lay = [0 1 3 10];   %probe layout
P = length(probe_lay);
delta_t = 1/omega/V; %virtual sampling interval

%signal generation
K_actual = 4;
k_initial = 10;

f = sort(randperm(V/2-1,K_actual)'*omega+randperm(omega/4,K_actual)'*3+normrnd(0,1,K_actual,1));%3 is not common divisor of rotating fre
f(3) = f(2)+1;

A = rand(length(f),1)*0.2+1;                  %a_k
A(1) = rand(1,1)*0.2+0.3;
A(4) = rand(1,1)*0.2+0.4;
phase = rand(length(f),1)*pi/2-pi/4;              %phase

Q_max = 130;
t_max = zeros(Q_max*P,1);
for i = 0:Q_max-1
    for j = 1:length(probe_lay)
        t_max(i*P+j) = (V*i+probe_lay(j))*delta_t; %the index of the actual sampling times
    end
end
x_pure_max = zeros(size(t_max));
for i = 1:length(f)
    x_pure_max=x_pure_max+A(i)*sin(2*pi*f(i)*t_max+phase(i)); %virtual displacement
end
power_sig = mean(x_pure_max.^2);

SNR = 10;
Qs = [10:10:100];
N_mc = 504;
Ks_estimate = zeros(length(Qs),6);
for i_Q = 1:length(Qs)
    waitbar(i_Q/length(Qs))
    Q = Qs(i_Q);                %number of revs
    M = Q*P;       %before padding (actual samples)
    
    var_noi = power_sig/10^(SNR/10);
    K_estimate = zeros(1,3);
    Probability = K_estimate;
    tic
    parfor n_mc = 1:N_mc
        seed = n_mc+(N_mc-1)*i_Q;
        rng(2022+seed)
        x_pure = x_pure_max(1:Q*P);
        t = t_max(1:Q*P);
        noise = normrnd(0,sqrt(var_noi),length(x_pure),1);
        x = x_pure+noise;
        %RELAX
        [Thitas,Var_s,f_fft] = RELAX_BTT_no_syn(x,t,probe_lay,V,k_initial,round(100000/V)*V);
        Thitas_cal = zeros(3*k_initial,k_initial);%FRE, AMP, PHASE
        for k = 1:k_initial
            Thitas_cal(3*k-2,1:k) = Thitas(3*k-2,1:k);
            Thitas_cal(3*k-1,1:k) = sqrt(Thitas(3*k-1,1:k).^2+Thitas(3*k,1:k).^2);
            Thitas_cal(3*k,1:k) = atan((Thitas(3*k,1:k)./Thitas(3*k-1,1:k)));
        end
%         [f_relax_sort,temp_index] = sort(Thitas_cal(3*k_initial-2,1:k_initial));
%         A_relax_sort = Thitas_cal(3*k_initial-1,temp_index);
%         Phase_relax_sort = Thitas_cal(3*k_initial,temp_index);

        %AICs and BICs
        AICs = zeros(size(Var_s));
        BICs = AICs; %from order=0
        AICs(1) = M*log(Var_s(1))+2*P;
        BICs(1) = M*log(Var_s(1))+P*log(M);
        for i = 2:k_initial+1
            BICs(i) = M*log(Var_s(i))+(5*(i-1)+P)*log(M);
            AICs(i) = M*log(Var_s(i))+(10*(i-1)+2*P);
        end
%         figure()
%         plot(0:k_initial,AICs)
%         xlabel("Model order")
%         ylabel("Value")
%         figure()
%         plot(0:k_initial,BICs)
%         xlabel("Model order")
%         ylabel("Value")
        %PAL
        PALs = zeros(k_initial,1);%1:K
        for i = 1:k_initial
            temp = M*log(Var_s(i+1)/M);
            r_n = M*log(Var_s(1)/Var_s(i));
            rho_n = M*log(Var_s(i)/Var_s(k_initial+1));
            PALs(i) = temp+(5*(i))*log(5*(k_initial))*log(r_n+1)/log(rho_n+1);
        end
        K_AIC = find(AICs==min(AICs(1:end)))-1;
        K_BIC = find(BICs==min(BICs(1:end)))-1;
        K_PAL = find(PALs==min(PALs));
        K_estimate = K_estimate+[K_AIC K_BIC K_PAL]/N_mc; 
        if K_AIC == K_actual
            Probability =Probability+[1 0 0];
        end
        if K_BIC == K_actual
            Probability =Probability+[0 1 0];
        end
        if K_PAL == K_actual
            Probability =Probability+[0 0 1];
        end
    end
    Ks_estimate(i_Q,1:3) = K_estimate;
    Ks_estimate(i_Q,4:6)= Probability/N_mc;
end
save('Ks_estimate_Q.mat','f','A','phase','Qs','Ks_estimate')





    
    