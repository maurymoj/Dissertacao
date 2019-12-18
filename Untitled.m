clc;clear;

% ENTRADAS DO PROGRAMA
%n = 1:365; % Dias do ano já convertidos para número total
n = 222;%[17 47 75 105 135 162 198 228 258 288 318 344];
lat = 43.0; % Latitude local [º]
altit = 1000; % Altitude local [m]

% Dados da superfície
beta_sup = 45;
gamma_sup = (lat<0)*180;

rho_g = 0.2;    % Refletividade do solo

% Radiação solar diária global
H = 0;
I = 0;

% Programa em si

% Ângulos diários da orientação do sol
[delta omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declinação solar
                                            % omega_s = ângulo de por do
                                            % sol

omega_s = real(omega_s);% ATENÇÃO !!!!!!!!!!!!!

% Radiação solar diária global no topo da atmosfera
H_o = Rad_Ext_Dia(lat, n, delta, omega_s);

% Radiação solar horária global no topo da atmosfera
omega_1 = floor(-max(omega_s)/15)*15 : 15 : floor(max(omega_s)/15)*15;
omega_2 =  ceil(-max(omega_s)/15)*15 : 15 : ceil(max(omega_s)/15)*15;
omega = (omega_1 + omega_2)/2;

% Radiação solar horária global no topo da atmosfera
I_o = Rad_Ext_Hor(lat, n, omega_1, omega_2, delta, omega_s);

if I == 0
    R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, delta);
    
    [I r_t] = Rad_I_de_H(H, lat, n, omega, delta, omega_s);% AVALIAR PRECISÃO DO CALCULO DE r_t e soma ao longo do dia.
end

% Índice de claridade
k_t = I./I_o;
k_t(isnan(k_t)) = 0;    % Retira os valores não válidos da matriz
k_t(isinf(k_t)) = 0;

% Estimativa das frações difusa e direta da radiação solar horária global
[I_b I_d] = frac_I(I, k_t);

% Cálculo dos ângulos de elevação solar e de azimute do sol
[theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

rho_r1 = 1.0;
alpha_1 = 15;
rho_r2 = rho_r1;
alpha_2 = 10;
[I_t_c I_t] = Rad_Inc_Conc(I, I_b, I_d, beta_sup, alpha_s, rho_r1, alpha_1, rho_r2, alpha_2, rho_g, R_b);