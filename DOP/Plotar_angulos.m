mes = 12;

figure('color',[1 1 1])
subplot(1,2,1)
plot(omega,gamma_sup(mes,:))
hold on
grid on
plot(omega,gamma_track(mes,:,1))
plot(omega,gamma_track(mes,:,2))
plot(omega,gamma_track(mes,:,3))
plot(omega,gamma_track(mes,:,4))
plot(omega,gamma_track(mes,:,5))
legend('fixo','Ideal','Relogio Biax','Proposto Biax','Relogio Uni','Proposto Uni')
title(cidades{cid})

subplot(1,2,2)
plot(omega,abs(lat)*ones(length(omega),1))
hold on
grid on
plot(omega,beta_track(mes,:,1))
plot(omega,beta_track(mes,:,2))
plot(omega,beta_track(mes,:,3))
plot(omega,beta_track(mes,:,4))
plot(omega,beta_track(mes,:,5))
legend('Fixo','Ideal','Relogio Biax','Proposto Biax','Relogio Uni','Proposto Uni')