% C�digo para estimativa da raz�o entre a �rea do painel e da �rea total
% necess�ria para prevenir sombreamento entre arrays
% Par�metros do problema

% Cidades e respectivos valores de latitude e altitude
cidades = {'Manaus (3.1�S)','Garanhuns (8.9�S)','Brasilia (15.8�S)','Belo Horizonte (19.8�S)','Campo Grande (20.4�S)','Joacaba (27.2�S)','Sao Gabriel (30.3�S)'};
lat_cid = [-3.10          -8.88             -15.77           -19.80                 -20.43               -27.17         -30.33];
altit_cid = [34.36        869.21            1115.25          937.53                 544.51               525.25         120.6];

n = data2dia(21,6); % Dia do solsticio de inverno no hemisf�rio sul
omega = -45; % 9 a.m.
theta = 0:80; % Varia��o do �ngulo de inclina��o do m�dulo

figure('color',[1 1 1]);
hold all;

for cid=1:length(cidades)
    lat = lat_cid(cid);

% �ngulos di�rios da orienta��o do sol
    [delta,omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declina��o solar
                                            % omega_s = �ngulo de por do
                                            % sol

% C�lculo dos �ngulos de zenite, eleva��o solar e de azimute do sol
    [theta_z,alpha_s,gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

% Raz�o entra as �reas total com espa�amento e �rea dos m�dulos de um array
    A_Rat = cosd(theta)+sind(theta)/tand(alpha_s)*cosd(180-gamma_s);
    
    plot(theta, A_Rat);
end

grid on;
legend(cidades,'Location','NorthWest');
xlabel('Tilt angle(�)')
ylabel('Area ratio A_r')
