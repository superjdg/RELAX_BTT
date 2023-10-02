%20221208
%@article{wang2023min, 
% title={Min-max Probe Placement and Extended Relaxation Estimation Method for Processing Blade Tip Timing Signals}, 
% author={Wang, Zengkun and Petre, Stoica and Dave, Zachariah and Prabhu, Babu and Zhibo, Yang}, 
% journal={IEEE TRANSACTIONS ON INSTRUMENTATION AND MEASUREMENT}, 
% year={2023}}
clear all
close all
rng(666)

f_r = 100;
n_vp = 24;
n_p = 4;
fs_v = f_r*n_vp;
N_rev=100;
t = [0:1/fs_v:(n_vp*N_rev-1)/fs_v]';
f = 0:fs_v/(n_vp*N_rev):(n_vp*N_rev-1)*fs_v/(n_vp*N_rev); %frequecy of FFT
fres = floor(rand(1,4)*fs_v/2);
a = [1 1 1 1];
phase = rand(1,4)*2*pi;
x = zeros(size(t));
for i = 1:length(fres)
    x = x+a(i)*sin(2*pi*fres(i)*t+phase(i));
end
layouts = [ 0 1 3 10; 0 2 5 15; 0 1 9 17;0 6 12 18;];
x_ffts = zeros(size(layouts,1),length(x)/2);
t_ffts = zeros(size(layouts,1),length(x)/2);

for idx_layout = 1:size(layouts,1)
    layout = layouts(idx_layout,:);
    t_index = zeros(size(t));
    for i = 0:N_rev-1
        t_index(layout+n_vp*i+1) = 1;
    end
    x_zero = t_index.*x;
    x_fft = abs(fft(x_zero));
    t_fft = abs(fft(t_index));
    x_fft = x_fft(1:length(x_fft)/2)/n_p/N_rev*2;
    t_fft = t_fft(1:length(t_fft)/2);
    x_ffts(idx_layout,:) = x_fft;
    t_ffts(idx_layout,:) = t_fft;
end

titles = ["(a) ","(b) ","(c) ","(d) "];
figure()
tiledlayout(size(layouts,1),1,'TileSpacing','tight')
for i = 1:size(layouts,1)
    nexttile
    h1 = plot(f(1:length(f)/2),x_ffts(i,:),'Color','#000000')
    hold on 
    for j = 1:length(fres)
        h2 = scatter(fres(j),x_fft(find(f==fres(j))),'x','black');
    end
    title(titles(i))
%     if i ==1
%         legend([h1,h2],{'${Y_{{\rm{pns}}}}$',['True frequency']}...
%    ,'NumColumns',1,'Interpreter','latex','FontName',"Times new roman","FontSize",6)
%         legend('boxoff')
%     end
    ylabel('Amplitude')
    if i ==4
        xlabel('Frequency/Hz')
    end
    ylim([0 1.2])
    set(gca,'ytick',[0 1])
end

%FFT of sampling sequence
titles = ["(a) ","(b) ","(c) ","(d) "];
figure()
tiledlayout(size(layouts,1),1,'TileSpacing','tight')
for i = 1:size(layouts,1)
    nexttile
    plot(f(1:length(f)/2),t_ffts(i,:),'Color','#000000',"Linewidth",2)
    hold on  
    title(titles(i))
    xlim([0 max(f(1:length(f)/2))+1])
    ylim([0 n_p*N_rev*1.2])

    ylabel('${|S|}$','Interpreter','latex')
    if i ==4
        xlabel('Frequency/Hz')
    end

end
