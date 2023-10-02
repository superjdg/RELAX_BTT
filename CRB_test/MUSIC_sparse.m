% 20221209
%@article{wang2023min, 
% title={Min-max Probe Placement and Extended Relaxation Estimation Method for Processing Blade Tip Timing Signals}, 
% author={Wang, Zengkun and Petre, Stoica and Dave, Zachariah and Prabhu, Babu and Zhibo, Yang}, 
% journal={IEEE TRANSACTIONS ON INSTRUMENTATION AND MEASUREMENT}, 
% year={2023}}
function [spectrum] = MUSIC_sparse(S,t,w, K)
%S snapshot matrix
%t time
%w frequency
%K number of components

    [U,D,V]=svd(S);
    signal_space = U(:,1:2*K);
%     signal_space = U(:,1:2*K+2);
    spectrum = zeros(size(w));
    for j = 1:length(w)%
        steer_vector = exp(1j*t*2*pi*w(j));
        steer_vector = steer_vector/norm(steer_vector);
        spectrum(j) = 1-abs((steer_vector'*signal_space)*(signal_space'*steer_vector));
    end
    spectrum = spectrum;
end


    
    
