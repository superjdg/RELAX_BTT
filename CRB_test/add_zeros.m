function [x,t] = add_zeros(t_zeros,x_nonzeros,t_nonzeros)
%@article{wang2023min, 
% title={Min-max Probe Placement and Extended Relaxation Estimation Method for Processing Blade Tip Timing Signals}, 
% author={Wang, Zengkun and Petre, Stoica and Dave, Zachariah and Prabhu, Babu and Zhibo, Yang}, 
% journal={IEEE TRANSACTIONS ON INSTRUMENTATION AND MEASUREMENT}, 
% year={2023}}
x = t_zeros;
t = x;
iter = 1;
for i = 1:length(t_zeros)
    if t_zeros(i)==1
        x(i) = x_nonzeros(iter);
        t(i) = t_nonzeros(iter);
        iter = iter+1;
    end
end
end