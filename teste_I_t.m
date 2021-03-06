% Script rascunho 2.0

clc;clear;

% VALORES PARA CADA CIDADE, NA ORDEM DE LATITUDE
lat_cid = [-3.10 -8.88 -15.77 -19.80 -20.43 -27.17 -30.33];
altit_cid = [39 841 1130 767 612 525 118];
H_cid = [16.82	16.32	16.68	17.60	16.85	15.08	16.74	19.87	20.23	18.78	20.45	17.14
         21.19	20.03	20.84	20.45	17.07	12.72	13.45	17.48	18.69	22.04	24.37	21.23
         21.94	20.71	19.76	20.04	17.69	14.33	17.37	21.15	23.26	24.44	20.84	18.50
         23.98	21.11	20.66	18.72	15.82	13.59	14.55	18.76	18.60	23.63	17.57	22.57
         22.75	21.30	19.60	20.86	16.52	14.69	14.13	18.19	19.89	23.11	25.95	23.59
         22.44	19.66	19.72	17.32	13.41	11.14	11.06	15.40	16.36	20.78	24.75	25.27
         23.81	21.28	19.53	14.41	12.04	9.23	9.84	14.53   17.54	18.26	24.59	25.60];

cid = 1; % Defini��o da cidade onde sera feita a simulacao

% ENTRADAS DO PROGRAMA
%n = 1:365; % Dias do ano j� convertidos para n�mero total
n = 75;
%n = [17 47 75 105 135 162 198 228 258 288 318 344]; % Dias m�dios de cada
%m�s
lat = 20.0; % Latitude local [�]
%lat = lat_cid(cid);
altit = 1000; % Altitude local [m] - somente necess�rio caso o modelo de
                % Perez et al for usado
%altit = altit_cid(cid);

% Dados da superf�cie
beta_sup = 0;
gamma_sup = (lat<0)*180;

rho_g = 0.2;    % Refletividade do solo

rastr = 0;  % O sistema possui sistema de rastreamento em dois eixos
conc = 1;   % O sistema possui sistema de concentra��o solar

% Radia��o solar di�ria global
H = 14.6323;

%H = [15.3290 20.0571 14.6323 15.3600 13.3548 13.4400 14.2839 15.6774 ...
%     16.8000 16.0258 17.2800 15.5613];

I = 0;
% Programa em si

% �ngulos di�rios da orienta��o do sol
[delta omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declina��o solar
                                            % omega_s = �ngulo de por do
                                            % sol

omega_s = real(omega_s);% ATEN��O !!!!!!!!!!!!!

% Radia��o solar di�ria global no topo da atmosfera
H_o = Rad_Ext_Dia(lat, n, delta, omega_s);

% Radia��o solar hor�ria global no topo da atmosfera
omega_1 = floor(-max(omega_s)/15)*15 : 15 : floor(max(omega_s)/15)*15;
omega_2 =  ceil(-max(omega_s)/15)*15 : 15 : ceil(max(omega_s)/15)*15;
omega = (omega_1 + omega_2)/2;

[theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);
% Radia��o solar hor�ria global no topo da atmosfera
I_o = Rad_Ext_Hor(lat, n, omega_1, omega_2, delta, omega_s);

%figure; hold on;

beta_sup = theta_z;
gamma_sup = gamma_s;

for i = 1:8

    if I == 0
        R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, delta);
    
        [I r_t] = Rad_I_de_H(H, lat, n, omega, delta, omega_s);% AVALIAR PRECIS�O DO CALCULO DE r_t e soma ao longo do dia.
    end

    % �ndice de claridade
    k_t = I./I_o;
    k_t(isnan(k_t)) = 0;    % Retira os valores n�o v�lidos da matriz
    k_t(isinf(k_t)) = 0;

    % Estimativa das fra��es difusa e direta da radia��o solar hor�ria global
    [I_b I_d] = frac_I(I, k_t);

    % C�lculo dos �ngulos de eleva��o solar e de azimute do sol


    I_t_perez = Rad_Inc(I, I_b, I_d, I_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup)

end