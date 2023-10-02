function [ x,residual] = Block_OMP( D, y, K ,e,f)
%@article{wang2023min, 
% title={Min-max Probe Placement and Extended Relaxation Estimation Method for Processing Blade Tip Timing Signals}, 
% author={Wang, Zengkun and Petre, Stoica and Dave, Zachariah and Prabhu, Babu and Zhibo, Yang}, 
% journal={IEEE TRANSACTIONS ON INSTRUMENTATION AND MEASUREMENT}, 
% year={2023}}
%=============================================
% This function realizes the Block orthogonal Matching Pursuit
% Input:
%       D : the dictionary you used (2*N columns for N frequencies:sine and cosine)
%       Y : the measured data,it must be the column vector
%       T : the number of parameters(to measure the sparsity,its 0-norm)
%       E : the tolerance of residual
% Output:
%       X : the sparse coefficient matrix
%       residual : the residual of the signal
% revised by Zengkun Wang [based on Zhibin Zhao's OMP code]
% zengkunwang@163.com
%=============================================

% Initialize the parameter matrix
N_f = size(D,2);
x = zeros(N_f,1);
N = length(y);
residual = y;
index = zeros(K,1);
W = zeros(N,2*K);
mean_var_old = mean(residual.^2);

for i = 1:K
    temp_norms = zeros(1,N_f/2);
    for j = 1:N_f/2
        phi_j = D(:,j*2-1:2*j);
        temp_norms(j) = residual'*phi_j*pinv(phi_j'*phi_j)*phi_j'*residual;
    end    
    index(i) = find(temp_norms==max(temp_norms));
    W(:,i*2-1:2*i) = D(:,index(i)*2-1:2*index(i));

    a_i = pinv(W'*W)*W'*y;
    residual = y - W * a_i; % Update the residual
    if mean(residual.^2)/mean_var_old < e % Stopping Rule
        break;
    end
    mean_var_old = mean(residual.^2);
end
for i = 1:length(index)
    x(2*index(i)-1:2*index(i)) = a_i(2*i-1:2*i);
end
end