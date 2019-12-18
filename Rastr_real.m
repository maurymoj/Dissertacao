mods_Rastr = {'fix_Int','fix_Int_Atr','lim_erro_incid'};
legend_Rastr = {'Algorithm based on d\omega','Algorithm based on delayed/ d\omega', 'Algorithm based on d\theta'};
%mods_Rastr = {'fix_int','lim_erro_beta','lim_erro_gamma','lim_erro_incid'};%,'prop_Aur','prop_Aur_2'};
%legend_Rastr = {'Algoritmo com base em d\omega','Algoritmo com base em d\beta','Algoritmo com base em d\gamma','Algoritmo com base em d\theta'};%,'Algoritmo ','prop Aur 2'};
crit_max = 80;  % Valor m�ximo para o crit�rio utilizado
crit_min = 0.1;   % Valor m�nimo para o crit�rio utilizado
crit_inc = 3;   % Valor do incremento para o crit�rio utilizado
C1 = 0.9:-0.2:0.1;
C2 = 0.1:0.2:0.9;
legend_Rastr_Aur = {'Algorithm based on d\omega','C1 = 0,9 e C2 = 0,1', 'C1 = 0,7 e C2 = 0,3', 'C1 = 0,5 e C2 = 0,5', 'C1 = 0,3 e C2 = 0,7', 'C1 = 0,1 e C2 = 0,9'};
d_omega = 0.1;    % Incremento nos valores do �ngulo da hora solar, [o].

cid = 7; % Defini��o da cidade onde sera feita a simulacao

% VALORES PARA CADA CIDADE, NA ORDEM DE LATITUDE - Somente s�o v�lidas
% latitudes com m�dulo menor ou igual a 66.5�
cidades = {'Manaus - AM','Garanhuns - PE','Brasilia - DF','Belo Horizonte - MG','Campo Grande - MS','Joacaba - SC','Sao Gabriel - RS'};
lat_cid = [-3.10          -8.88             -15.77           -19.80                 -20.43               -27.17         -30.33];
altit_cid = [34.36        869.21            1115.25          937.53                 544.51               525.25         120.6];

% ENTRADAS DO PROGRAMA

% CONSTANTES

rho_g = 0.2;    % Refletividade do solo

G_sc = 1367; % Constante solar

% Dias do ano j� convertidos para n�mero total
n = [17 47 75 105 135 162 198 228 258 288 318 344]; % Dias m�dios de cada m�s
n_dias_mes = [31 28 31 30 31 30 31 31 30 31 30 31]; % N�mero de dias do respectivo m�s


% Latitude local [�] 
lat = lat_cid(cid);

% Altitude local [m]
altit = altit_cid(cid);

% Radia��o solar di�ria global

H = H_cid(cid,:);

% GEOMETRIA DA RADIA��O SOLAR

% �ngulos di�rios da orienta��o do sol
[delta omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declina��o solar
                                            % omega_s = �ngulo de por do
                                            % sol
mods_Rastr = {'fix_Int','fix_Int_Atr','lim_erro_incid'};
legend_Rastr = {'Algorithm based on d\omega','Algorithm based on delayed/ d\omega', 'Algorithm based on d\theta'};
%mods_Rastr = {'fix_int','lim_erro_beta','lim_erro_gamma','lim_erro_incid'};%,'prop_Aur','prop_Aur_2'};
%legend_Rastr = {'Algoritmo com base em d\omega','Algoritmo com base em d\beta','Algoritmo com base em d\gamma','Algoritmo com base em d\theta'};%,'Algoritmo ','prop Aur 2'};
crit_max = 80;  % Valor m�ximo para o crit�rio utilizado
crit_min = 0.1;   % Valor m�nimo para o crit�rio utilizado
crit_inc = 3;   % Valor do incremento para o crit�rio utilizado
C1 = 0.9:-0.2:0.1;
C2 = 0.1:0.2:0.9;
legend_Rastr_Aur = {'Algorithm based on d\omega','C1 = 0.9 and C2 = 0.1', 'C1 = 0.7 and C2 = 0.3', 'C1 = 0.5 e C2 = 0.5', 'C1 = 0.3 and C2 = 0.7', 'C1 = 0.1 and C2 = 0.9'};
d_omega = 0.1;    % Incremento nos valores do �ngulo da hora solar, [o].