close all
%K 
figure()
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile
load Ks_estimate_snr
Ks_estimate_snr = Ks_estimate;
colors = ["#049CD8","#F2921D","#43B047","#E52521","#000000"];
marks = ['--s',"--+","--^"];
hold on
for j = 1:3
    plot(SNRs,Ks_estimate_snr(:,j),marks(j),'Color',colors(j));
end
box on
ylabel("mean$(\hat K)$",'Interpreter','latex')
xlabel("SNR/dB")
legend(["AIC","BIC","PAL"],'Location','northwest')
legend('boxoff')
ylim([0 10])
title("(a)")

nexttile
load Ks_estimate_Q
Ks_estimate_Q = Ks_estimate;
hold on
for j = 1:3
    plot(Qs,Ks_estimate_Q(:,j),marks(j),'Color',colors(j));
end
box on
xlabel("$Q$",'Interpreter','latex')
ylim([0 10])
xlim([10 100])
title("(b)")

%probability
figure()
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
for j = 1:3
    plot(SNRs,Ks_estimate_snr(:,j+3),marks(j),'Color',colors(j));
end
box on
ylabel("Probability",'Interpreter','latex')
xlabel("SNR/dB")
legend(["AIC","BIC","PAL"],'Location','northwest')
legend('boxoff')
ylim([0 1])

nexttile
title("(b)")
hold on
for j = 1:3
    plot(Qs,Ks_estimate_Q(:,j+3),marks(j),'Color',colors(j));
end
box on
xlabel("$Q$",'Interpreter','latex')
ylim([0 1])
xlim([10 100])

% set(gca,'Ytick',[0,10])
