clear all
close all
%% vs SNR
load MSEs_snr_white
RMSEs = sqrt(MSEs);
RMSE_F_SNR = zeros(length(SNRs),4);
for i = 1
%     nexttile
    hold on
    for j = 1:4
        temp1 = RMSEs(:,(j-1)*3+i);
        temp2 = zeros(size(temp1,1)/4,1);
        for k = 1:4
            temp2 = temp2+temp1(k:4:end)/4/f(k);
        end
        RMSE_F_SNR(:,j) = temp2;
    end
end

colors = ["#049CD8","#F2921D","#43B047","#E52521","#000000"];
marks = ['-s',"-+","-^","-diamond"];
hold on
for j = 1:4
    plot(SNRs,RMSE_F_SNR(:,j),marks(j),'Color',colors(j));
end
box on
ylabel("Normalized RMSE")
xlabel("SNR/dB")
ylim([min(min(RMSE_F_SNR)),max(max(RMSE_F_SNR))]);
set(gca,'YScale','log')
legend(["RELAX","Block-OMP","MUSIC","CRB"],'Location','northwest')
legend('boxoff')
title('(a)')


%% VS Q
% clear all
load MSEs_Q_white
RMSEs = sqrt(MSEs);
RMSE_F_Q = zeros(length(Qs),4);
for i = 1
%     nexttile
    hold on
    for j = 1:4
        temp1 = RMSEs(:,(j-1)*3+i);
        temp2 = zeros(size(temp1,1)/4,1);
        for k = 1:4
            temp2 = temp2+temp1(k:4:end)/4/f(k);
        end
        RMSE_F_Q(:,j) = temp2;
    end
end


colors = ["#049CD8","#F2921D","#43B047","#E52521","#000000"];
marks = ['-s',"-+","-^","-diamond"];
figure()
hold on
for j = 1:4
    plot(Qs,RMSE_F_Q(:,j),marks(j),'Color',colors(j));
end
box on
ylabel("Normalized RMSE")
xlabel("$Q$",'Interpreter','latex')
ylim([min(min(RMSE_F_Q)),max(max(RMSE_F_Q))]);
xlim([10 130])
set(gca,'YScale','log')
legend(["RELAX","Block-OMP","MUSIC","CRB"],'Location','northwest')
legend('boxoff')
title('(b)')
