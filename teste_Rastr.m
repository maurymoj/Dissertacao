% Script rascunho 2.5

clc;clear;

% VALORES PARA CADA CIDADE, NA ORDEM DE LATITUDE
lat_cid = [-3.10 -8.88 -15.77 -19.80 -20.43 -27.17 -30.33];
altit_cid = [39 841 1130 767 612 525 118];
% H_cid = [16.82	16.32	16.68	17.60	16.85	15.08	16.74	19.87	20.23	18.78	20.45	17.14
%          21.19	20.03	20.84	20.45	17.07	12.72	13.45	17.48	18.69	22.04	24.37	21.23
%          21.94	20.71	19.76	20.04	17.69	14.33	17.37	21.15	23.26	24.44	20.84	18.50
%          23.98	21.11	20.66	18.72	15.82	13.59	14.55	18.76	18.60	23.63	17.57	22.57
%          22.75	21.30	19.60	20.86	16.52	14.69	14.13	18.19	19.89	23.11	25.95	23.59
%          22.44	19.66	19.72	17.32	13.41	11.14	11.06	15.40	16.36	20.78	24.75	25.27
%          23.81	21.28	19.53	14.41	12.04	9.23	9.84	14.53   17.54	18.26	24.59	25.60];

% ENTRADAS DO PROGRAMA

rad = 3; % Definicao da fonte dos dados de radia��o solar:
         %              1 = Radia��o solar hor�ria
         %              2 = Radia��o solar di�ria
         %              3 = Radia��o solar estimada para dia claro
         
cid = 1; % Defini��o da cidade onde sera feita a simulacao

%n = 1:365; % Dias do ano j� convertidos para n�mero total
%n = 75;
n = [17 47 75 105 135 162 198 228 258 288 318 344]; % Dias m�dios de cada m�s
%lat = 50.0; % Latitude local [�]
lat = lat_cid(cid);
altit = 612; % Altitude local [m] - somente necess�rio caso o modelo de
                % Perez et al for usado
%altit = altit_cid(cid);

% Dados da superf�cie
beta_sup = abs(lat);
gamma_sup = (lat<0)*180;

rho_g = 0.2;    % Refletividade do solo

rastr = 0;  % O sistema possui sistema de rastreamento em dois eixos
conc = 0;   % O sistema possui sistema de concentra��o solar

% Radia��o solar di�ria global

% H = 14.6323;
H = [15.3290 20.0571 14.6323 15.3600 13.3548 13.4400 14.2839 15.6774 ...
     16.8000 16.0258 17.2800 15.5613];

I = 0;
% Programa em si

% GEOMETRIA DA RADIA��O SOLAR

% �ngulos di�rios da orienta��o do sol
[delta omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declina��o solar
                                            % omega_s = �ngulo de por do
                                            % sol

omega_s = real(omega_s);% ATEN��O !!!!!!!!!!!!!

omega_1 = floor(-max(omega_s)/15)*15 : 15 : floor(max(omega_s)/15)*15;
omega_2 =  ceil(-max(omega_s)/15)*15 : 15 : ceil(max(omega_s)/15)*15;

% RADIA��O SOLAR NO TOPO DA ATMOSFERA

% Radia��o solar di�ria global no topo da atmosfera
H_o = Rad_Ext_Dia(lat, n, delta, omega_s);

% Radia��o solar hor�ria global no topo da atmosfera
I_o = Rad_Ext_Hor(lat, n, omega_1, omega_2, delta, omega_s);

% C�LCULO DA RADIA��O SOLAR GLOBAL, DIRETA E DIFUSA EM SUPERF�CIE
% HORIZONTAL

if rad == 2
    omega = (omega_1 + omega_2)/2;
    
    % C�lculo dos �ngulos de eleva��o solar e de azimute do sol
    [theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);
    
    % C�lculo da radia��o solar global hor�ria a partir da radia��o solar global di�ria
    [I r_t] = Rad_I_de_H(H, lat, n, omega, delta, omega_s);% AVALIAR PRECIS�O DO CALCULO DE r_t e soma ao longo do dia.

    % �ndice de claridade
    k_t = I./I_o;
    k_t(isnan(k_t)) = 0;    % Retira os valores n�o v�lidos da matriz
    k_t(isinf(k_t)) = 0;

    % Estimativa das fra��es difusa e direta da radia��o solar hor�ria global
    [I_b I_d] = frac_I(I, k_t);
elseif rad == 3
    omega = floor(-max(omega_s)/15)*15 : 1 : ceil(max(omega_s)/15)*15;
    
    % C�lculo dos �ngulos de eleva��o solar e de azimute do sol
    [theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);  
    
    clima = 1;% Definicao do clima da regi�o estudada:
              %              1 - Clima tropical
              %              2 - Verao de latitude media
              %              3 - Verao sub artico
              %              4 - Inverno de latitude media
    
    B = 2*pi/360*(n-1)*360/365;	% Variavel auxiliar
    
    % Radiacao extraterrestre em superficie perpendicular a radiacao,
    % conforme Duffie e Beckman (2006) e Iqbal (1983)

	G_on = G_sc*(1.000110 + 0.034221*cos(B) + 0.001280*sin(B)...	
		   + 0.000719*cos(2*B) + 0.000077*sin(2*B));
       
    [G_h G_b G_d] = Rad_Dia_Claro(theta_z, altit, clima, G_on);
    
end

% C�LCULO DA RADIA��O EM SUPERF�CIE COM ORIENTA��O ARBITR�RIA

% Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
% sobre uma superf�cie horizontal
R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, delta);
R_b(R_b<0)=0; %INTEGRACAO DE Rb COM HORA CONTENDO POR DO SOL NASCER DO SOL

% Radia��o em superf�cie arbitr�ria de acordo com o modelo de Perez et al
I_t_perez = Rad_Inc(I, I_b, I_d, I_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);

if ~rastr && conc
    rho_r1 = 0.7;       % Refletividade dos espelhos
    alpha_1 = 73;
    rho_r2 = rho_r1;
    alpha_2 = 17;
    [I_t_c I_t] = Rad_Inc_Conc(I, I_b, I_d, beta_sup, alpha_s, rho_r1, alpha_1, rho_r2, alpha_2, rho_g, R_b);
    g_c = sum(I_t_c,2)./sum(I_t,2)
    [sum(sum(I_t_c,2)) sum(sum(I_t,2))]
elseif rastr && ~conc
    % Inserir trecho de altera��o/inser��o do crit�rio de reorienta��o ?
    [I_t_r beta_sup gamma_sup] = Rad_Inc_Rast(I, I_b, I_d, I_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, 'ideal', 'inc_Dia', 10);
    sum(sum(I_t_r,2))./sum(sum(I_t_perez,2))
end