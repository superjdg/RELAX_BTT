% 20221208
%@article{wang2023min, 
% title={Min-max Probe Placement and Extended Relaxation Estimation Method for Processing Blade Tip Timing Signals}, 
% author={Wang, Zengkun and Petre, Stoica and Dave, Zachariah and Prabhu, Babu and Zhibo, Yang}, 
% journal={IEEE TRANSACTIONS ON INSTRUMENTATION AND MEASUREMENT}, 
% year={2023}}
clear all
close all

n_vp = 24;
n_p = 4;
layouts = combntns([0:n_vp-1],n_p);
layouts = layouts-layouts(:,1);%location of first probe is 0
uni_layouts = unique(layouts,'rows',"stable");

max_norms = zeros(size(uni_layouts,1),1);
for g = 1:size(uni_layouts,1)%
    waitbar(g/size(uni_layouts,1))
    layout = uni_layouts(g,:);
    cks = zeros(n_vp,1);
    for i = 1:n_vp
        z = zeros(n_p,1);
        for j = 1:n_p
            z(j) = exp(-2*pi*1j*(i-1)/n_vp*(layout(j)));
        end
        cks(i) = sum(z);
    end
    max_norms(g) = max(abs(cks(2:end))/abs(cks(1)));
end

[select_index] = find(max_norms-min(max_norms)<0.0000001/min(max_norms));
layouts_selec = uni_layouts(select_index,:);

intervals = zeros(size(layouts_selec,1),n_p);
for i = 1:size(layouts_selec,1)
    temp_layout = layouts_selec(i,:);
    for j = 1:1
        intervals(i,(j-1)*n_p+1:j*n_p) = [temp_layout(j+1:end)-temp_layout(1:n_p-j)...
            temp_layout(1:j)+n_vp-temp_layout(n_p-j+1:end)];
    end
end
layouts_selec_unique = layouts_selec;
for i = 1:size(layouts_selec_unique,1)
    temp_interval = intervals(i,:);
    if sum(temp_interval) == 0
        continue
    end
    for j = 2:n_p
        temp_interval_circular = [temp_interval(j:end) temp_interval(1:j-1)];
        for k = i:size(layouts_selec_unique,1)
            if sum(abs(intervals(k,:)-temp_interval_circular))==0 ...
                    ||sum(abs(fliplr(intervals(k,:))-temp_interval_circular)) ==0
                layouts_selec_unique(k,:) = zeros(1,n_p);
                continue
            end
        end
    end
end
layouts_selec_unique(all(layouts_selec_unique==0,2),:) = [];
layouts_selec_unique

[worst_index] = max(max_norms);
worst_layout = uni_layouts(worst_index,:);

N_rev = 50;
for i = 1:size(layouts_selec_unique,1)
    temp_layout = layouts_selec_unique(i,:);
    sample_sequence = zeros(1,N_rev*n_vp);
    for j = 1:N_rev
        sample_sequence((j-1)*n_vp+1+temp_layout)=1;
    end
    padded_zeros = zeros(1,9973-N_rev*n_vp);
    sample_sequence_padded = [sample_sequence padded_zeros];
end


