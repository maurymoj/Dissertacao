% Código para estimativa da razão entre a área do painel e da área total
% necessária para prevenir sombreamento entre arrays
% Parâmetros do problema

% Cidades e respectivos valores de latitude e altitude
cidades = {'Manaus (3.1°S)','Garanhuns (8.9°S)','Brasilia (15.8°S)','Belo Horizonte (19.8°S)','Campo Grande (20.4°S)','Joacaba (27.2°S)','Sao Gabriel (30.3°S)'};
lat_cid = [-3.10          -8.88             -15.77           -19.80                 -20.43               -27.17         -30.33];
altit_cid = [34.36        869.21            1115.25          937.53                 544.51               525.25         120.6];

n = data2dia(21,6); % Dia do solsticio de inverno no hemisfério sul
omega = -45; % 9 a.m.
theta = 0:80; % Variação do ângulo de inclinação do módulo

figure('color',[1 1 1]);
hold all;

for cid=1:length(cidades)
    lat = lat_cid(cid);

% Ângulos diários da orientação do sol
    [delta,omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declinação solar
                                            % omega_s = ângulo de por do
                                            % sol

% Cálculo dos ângulos de zenite, elevação solar e de azimute do sol
    [theta_z,alpha_s,gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

% Razão entra as áreas total com espaçamento e área dos módulos de um array
    A_Rat = cosd(theta)+sind(theta)/tand(alpha_s)*cosd(180-gamma_s);
    
    plot(theta, A_Rat);
end

grid on;
legend(cidades,'Location','NorthWest');
xlabel('Tilt angle(°)')
ylabel('Area ratio A_r')
