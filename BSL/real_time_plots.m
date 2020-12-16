function void = real_time_plots(theta_true, thetas, j, accept, k, prior,prop,likelihood,bestlikelihood)

ttt = tiledlayout(4,1);
title(ttt,['Number of steps: ',num2str(j),'  Acceptance rate: ',num2str(accept/j), ' Current likelihood: ', num2str(likelihood), ' Best likelihood: ', num2str(bestlikelihood),' Probability to accept: ', num2str(exp(likelihood-bestlikelihood)) ])
nexttile
hold on
plot(thetas(1,1:j),'o')
plot(j,prop(1),'r*','Markersize',24)
ylim([prior(1,1) prior(1,2)])
yline(theta_true(1))
ylabel("Magnitude")
title("T")

nexttile
hold on
plot(pow2db(thetas(2,1:j)),'o')
plot(j,pow2db(prop(2)),'r*','Markersize',24)
ylim([pow2db(prior(2,1)) pow2db(prior(2,2))])
yline(pow2db(theta_true(2)))
ylabel("dB")
title("G0")

nexttile
hold on
plot(thetas(3,1:j),'o')
plot(j,prop(3),'r*','Markersize',24)
ylim([prior(3,1) prior(3,2)])
yline(theta_true(3))
ylabel("Magnitude")
title("\lambda")

nexttile
hold on
plot(thetas(4,1:j),'o')
plot(j,prop(4),'r*','Markersize',24)
ylim([prior(4,1) prior(4,2)])
yline(theta_true(4))
ylabel("Magnitude")
title("\sigma")

pause(.00000001);
if j ~=k-1
    delete(ttt)
end
void = 0;
end