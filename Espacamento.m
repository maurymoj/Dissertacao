lat=3.1;
alpha_s=57.24;
beta=1:30;
L=2.6;

lat = [3.1 8.88 15.77 19.80 20.43 27.17 30.33 ];
alpha_s = [62.85 57.24 50.5050 46.54 45.9215 39.2722 36.1489];

D1 = L*cosd(beta)+L*sind(beta)/tand(alpha_s(1));
D2 = L*cosd(beta)+L*sind(beta)/tand(alpha_s(2));
D3 = L*cosd(beta)+L*sind(beta)/tand(alpha_s(3));
D4 = L*cosd(beta)+L*sind(beta)/tand(alpha_s(4));
D5 = L*cosd(beta)+L*sind(beta)/tand(alpha_s(5));
D6 = L*cosd(beta)+L*sind(beta)/tand(alpha_s(6));
D7 = L*cosd(beta)+L*sind(beta)/tand(alpha_s(7));

figure
plot(beta,D1)
hold all
plot(beta,D2)
plot(beta,D3)
plot(beta,D4)
plot(beta,D5)
plot(beta,D6)
plot(beta,D7)
legend(lat(1),lat(2),lat(3),lat(4),lat(5),lat(6),lat(7))