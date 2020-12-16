thetas_plot = thetas(:,1:j);

tt = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');

len = length(thetas_plot(1,1:j));

linesize = 2;
fontsize = 8;
truecolor = '#32CD32';

burnin = 200;

T_true = theta_true(1);
T_start = thetas(1,1);

G0_true = theta_true(2);
G0_start = thetas(2,1);

lambda_true = theta_true(3);
lambda_start = thetas(3,1);

sigma_N_true = theta_true(4);
sigma_start = thetas(4,1);

nexttile

[f, x] = ksdensity(thetas_plot(1,burnin:len/blocks*(kk)),linspace(prior(1,1), prior(1,2),100000));
area(x*1e9,f,'FaceColor','#bbbbbb','LineStyle','None')
MMSE_T = mean(thetas(1,burnin:j));
xline(MMSE_T*1e9,'r','LineWidth',linesize)
xline(T_true*1e9,'--','Color',truecolor,'LineWidth',linesize)
subtitle('T \times 10^{-9}')
xlim([prior(1,1)*1e9 prior(1,2)*1e9])
%     xticks([prior(1,1) 4.5e-9 8e-9 11.5e-9 prior(1,2)])
%set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])


legend( "Approx. posterior",'MMSE',"True value")

thetas_plot(2,:) = pow2db(thetas_plot(2,:));

nexttile
[f, x] = ksdensity(thetas_plot(2,burnin:len/blocks*(kk)),pow2db(linspace(prior(2,1), prior(2,2),100000)));
area(x,f,'FaceColor','#bbbbbb','LineStyle','None')
MMSE_G0 = mean(thetas_plot(2,burnin:j-1));
xline(MMSE_G0,'r','LineWidth',linesize)
xline(pow2db(G0_true),'--','Color',truecolor,'LineWidth',linesize)
subtitle('G_0 [dB]')
xlim([pow2db(prior(2,1)) pow2db(prior(2,2))])
%set(gca,'xticklabel',[])
xticks([-90 -85 -80 -75])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

nexttile

[f, x] = ksdensity(thetas_plot(3,burnin:len/blocks*(kk)),linspace(prior(3,1), prior(3,2),100000));
area(x*1e-9,f,'FaceColor','#bbbbbb','LineStyle','None')
MMSE_lambda = mean(thetas(3,burnin:j));
xline(MMSE_lambda*1e-9,'r','LineWidth',linesize)
xline(lambda_true*1e-9,'--','Color',truecolor,'LineWidth',linesize)
subtitle('\lambda \times 10^{9}')
xlim([5e6*1e-9 4e9*1e-9])
%set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])


nexttile

[f, x] = ksdensity(thetas_plot(4,burnin:len/blocks*(kk)),linspace(prior(4,1), prior(4,2),100000));
area(x.^2*1e9,f,'FaceColor','#bbbbbb','LineStyle','None')
MMSE_sigma_N = mean(thetas(4,burnin:j));
xline(MMSE_sigma_N^2*1e9,'r','LineWidth',linesize)
xline(sigma_N_true^2*1e9,'--','Color',truecolor,'LineWidth',linesize)
subtitle('\sigma_N^2 \times 10^{-9}')
xlim([prior(4,1)^2*1e9 prior(4,2)^2*1e9])
%set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

