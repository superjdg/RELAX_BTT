function [Thitas,Var_s,f] = RELAX_BTT_no_syn(x_no_zeros,t_no_zeros,probe_layout,n_vp,K,N_paded)
%RELAX 
%@article{wang2023min, 
% title={Min-max Probe Placement and Extended Relaxation Estimation Method for Processing Blade Tip Timing Signals}, 
% author={Wang, Zengkun and Petre, Stoica and Dave, Zachariah and Prabhu, Babu and Zhibo, Yang}, 
% journal={IEEE TRANSACTIONS ON INSTRUMENTATION AND MEASUREMENT}, 
% year={2023}}
n_p = length(probe_layout);%probe Layout 
n_rev = length(x_no_zeros)/n_p;
fr = 1/mean(t_no_zeros(1+n_p:n_p:end)-t_no_zeros(1:n_p:end-n_p));
fs = n_vp*fr;
N = n_vp*n_rev;
N_a = n_p*n_rev;
% N_paded = 1200*N;
zeros_added = zeros(N_paded-N,1);

t_zero = zeros(N,1);
for i = 1:n_rev
    t_zero((i-1)*n_vp+probe_layout+1) = 1;
end

% BICs = zeros(K,1); %BIC values
Thitas = zeros(3*K,K); %for the frequencies and amplitudes of sine and cosine
Var_s = zeros(K+1,1);


f = 0:fs/N_paded:(N_paded-1)*fs/N_paded; %frequecy of FFT


% k = 0 %initialization
variance = x_no_zeros'*x_no_zeros;
Var_s(1) = variance;
residual = x_no_zeros;
x_k = x_no_zeros; %initialization of residual


Phi = [];
A_2k = [];
F = [];
for k = 1:K
    Phi = [Phi zeros(N_a,2)];%add_one_more_frequency
    A_2k = [A_2k zeros(2,1)];
    F = [F 0];
    for iter = 1:200
        temp_index = [k 1:k-1];%0 denote the DC components
        variance_old = variance;
        for jj = 1:length(temp_index)
            kk = temp_index(jj);
%             disp([k,iter,kk])
            
            [x_k_zeros,~] = add_zeros(t_zero,x_k,t_no_zeros);
            residual_zeros_tail = [x_k_zeros;zeros_added];
            x_fft_temp = fft(residual_zeros_tail);
%             figure()
%             plot(f,abs(x_fft_temp))
            x_fft_temp_half = abs(x_fft_temp(1:N_paded/2));
            
            f_index = find(x_fft_temp_half==max(x_fft_temp_half));
            f_selected = f(f_index);
            phi_k = [sin(2*pi*f_selected*t_no_zeros) cos(2*pi*f_selected*t_no_zeros) ];
            a_2k = pinv(phi_k'*phi_k)*phi_k'*x_k;
            F(kk) = f_selected;
            Phi(:,kk*2-1:2*kk) = phi_k;
            A_2k(:,kk) = a_2k;
            next_kk = mod(kk,k)+1;%the one to predict in next step
            x_k = x_no_zeros-[Phi(:,1:(next_kk-1)*2) Phi(:,next_kk*2+1:end)]...
                *reshape([A_2k(:,1:next_kk-1) A_2k(:,next_kk+1:end)],[2*(k-1),1]);        
        end
        residual = x_no_zeros-Phi*reshape(A_2k,[2*k,1]); %for compute converge
        variance = residual'*residual;%
%         variance
%         variance_old
        if abs((variance-variance_old)/variance_old)<0.0001%
            Var_s(k+1) = variance;
            Thitas(k*3-2:k*3,1:k) = [F; A_2k];
            break
        end
        Var_s(k+1) = variance;
        Thitas(k*3-2:k*3,1:k) = [F; A_2k];
    end
    x_k=residual;
end

end

