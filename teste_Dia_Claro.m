% Script rascunho 2.5

clc;clear;

G_sc = 1367; %[W/m^2]

% VALORES PARA CADA CIDADE, NA ORDEM DE LATITUDE
lat_cid = [-3.10 -8.88 -15.77 -19.80 -20.43 -27.17 -30.33];
altit_cid = [39 841 1130 767 612 525 118];

cid = 1; % Defini��o da cidade onde sera feita a simulacao

% ENTRADAS DO PROGRAMA

n = [17 47 75 105 135 162 198 228 258 288 318 344]; % Dias m�dios de cada m�s

lat = lat_cid(cid); % Latitude local [�]
altit = altit_cid(cid); % Altitude local [m] - somente necess�rio caso o 
                        % modelo de Perez et al seja usado

% Dados da superf�cie
beta_sup = abs(lat);
gamma_sup = (lat<0)*180;

rho_g = 0.2;    % Refletividade do solo

rastr = 0;  % O sistema possui sistema de rastreamento em dois eixos
conc = 0;   % O sistema possui sistema de concentra��o solar


% Programa em si


% GEOMETRIA DA RADIA��O SOLAR

% �ngulos di�rios da orienta��o do sol
[delta omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declina��o solar
                                            % omega_s = �ngulo de por do
                                            % sol

omega_s = real(omega_s);% ATEN��O !!!!!!!!!!!!!

omega_1 = floor(-max(omega_s)/15)*15 : 15 : floor(max(omega_s)/15)*15;
omega_2 =  ceil(-max(omega_s)/15)*15 : 15 : ceil(max(omega_s)/15)*15;
%omega = (omega_1 + omega_2)/2;

omega = floor(-max(omega_s)/15)*15 : 1 : ceil(max(omega_s)/15)*15;

% C�lculo dos �ngulos de eleva��o solar e de azimute do sol
[theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

% RADIA��O SOLAR NO TOPO DA ATMOSFERA

% Radia��o solar di�ria global no topo da atmosfera
H_o = Rad_Ext_Dia(lat, n, delta, omega_s);

% Radia��o solar hor�ria global no topo da atmosfera
I_o = Rad_Ext_Hor(lat, n, omega_1, omega_2, delta, omega_s);

B = 2*pi/360*(n-1)*360/365;	% Variavel auxiliar
    
% Calculo radiacao de dia claro para superficie horizontal

G_on = G_sc*(1.000110 + 0.034221*cos(B) + 0.001280*sin(B)...	% Radiacao extraterrestre em superficie perpendicular a radiacao, conforme Duffie e Beckman (2006) e Iqbal (1983)
		+ 0.000719*cos(2*B) + 0.000077*sin(2*B));

[G_t G_b G_d] = Rad_Dia_Claro(theta_z, A, clima, G_on);

% C�LCULO DA RADIA��O EM SUPERF�CIE COM ORIENTA��O ARBITR�RIA

% Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
% sobre uma superf�cie horizontal
R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, delta);
R_b(R_b<0)=0; %INTEGRACAO DE Rb COM HORA CONTENDO POR DO SOL NASCER DO SOL

% Radia��o em superf�cie arbitr�ria de acordo com o modelo de Perez et al
I_t_perez = Rad_Inc(I, I_b, I_d, I_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);