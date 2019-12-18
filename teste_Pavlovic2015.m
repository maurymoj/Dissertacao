clear;clc;

G_sc = 1366;

n = 47;%[17 47 75 105 135 162 198 228 258 288 318 344];
lat = 40.0; % Latitude local [º]
altit = 0.194;% Altitude local [km]

% Dados da superfície
beta_sup = 50;
gamma_sup = (lat<0)*180;

rho_g = 0.2;    % Refletividade do solo

% Ângulos diários da orientação do sol
[delta omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declinação solar
                                            % omega_s = ângulo de por do
                                            % sol

omega_s = real(omega_s);% ATENÇÃO !!!!!!!!!!!!!

% Radiação solar horária global no topo da atmosfera
omega_1 = ((5-12)*15 : 15 : (18-12)*15);
omega_2 = ((6-12)*15 : 15 : (19-12)*15);
omega = (omega_1 + omega_2)/2;
%omega = ((6:19)-12)*15+7.5;

[theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

G_o = G_sc*(1 + 0.033*cosd(360*n/365)).*cosd(theta_z);
%G   = 0.9*G_o;
G   = [25 150 325 500 675 800 850 875  825 750  600 425 225 0];
%           6   7   8   9  10  11  12  13  14   15  16   17  18 19
k_t = G./G_o;
k_t(isnan(k_t)) = 0;    % Retira os valores não válidos da matriz
k_t(isinf(k_t)) = 0;

[G_b G_d] = frac_I(G, k_t);

%G_b = G_sc*( (1-a*altit)*C_1.^(AM.^(C_2)) + a*altit ); % Para esta equação é necessário o valor da altitude em kilometros !

R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, delta);
R_b(omega_1 < -omega_s) = 0;
R_b(omega_2 > omega_s)  = 0;
R_b(R_b < 0) = 0;

G_dir_col = G_b.*(sind(alpha_s + beta_sup)./cosd(theta_z));

G_dif_col = G_d.*(1 + cosd(beta_sup))/2;

% Radiação refletida pelo solo
G_ref_gr = rho_g*G.*(1-cosd(beta_sup))/2;

alpha_1 = 38;
alpha_2 = 5;
rho_r1 = 1.0;
rho_r2 = rho_r1;

% Radiação refletida pelo refletor 1
Qui = beta_sup + 2*alpha_1 - alpha_s;

G_ref_r1 = max(rho_r1*(G_b./cosd(theta_z)).*sind(Qui).*sind(alpha_s - alpha_1),0);

% Radiação refletida pelo refletor 
Tau = alpha_s + 2*alpha_2 - beta_sup;

G_ref_r2 = max(rho_r2*(G_b./cosd(theta_z)).*sind(Tau).*cosd(alpha_s + alpha_2),0);

% Radiação total incidente sobre o módulo
G_t_conc = G_dir_col + G_dif_col + G_ref_gr + G_ref_r1 + G_ref_r2;

G_t = G_dir_col + G_dif_col + G_ref_gr;

G_t(G_t < 0) = 0;
G_t(alpha_s < 0) = 0;

G_t_conc(G_t_conc < 0) = 0;
G_t_conc(alpha_s < 0) = 0;

figure;
plot(omega, G_t, omega, G_t_conc,'r');
grid on;