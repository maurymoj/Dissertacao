%% Etapa 1 - Padronizacao dos dados medidos na estacao meteorologica de Vicosa

% DADOS DA CIDADE DE VI�OSA - Conforme cadastro de localidades do IBGE 
% divulgado em 2010
lat = -20.75;       % [�]
%alt = 691.61;      % [m]
long = 42.88;       % [�] longitude medida em graus a oeste (de 0� a 360�)
long_std = 47.89;   % [�]

% Carregando dados de radia��o para Vi�osa
Dados = load('DadosReaisVic.txt');
                                   % Cada linha desta matriz cont�m os
                                   % dados para um determinado momento de
                                   % um dia em particular, as colunas est�o
                                   % organizadas da seguinte forma: 
                                   % [dia, m�s, ano, hora, minuto, radia��o 
                                   % m�dia em W/m�, e a ultima coluna tem a
                                   % radi�ncia para o intervalo em mJ/m�.
                                   
                                   % Os valores incompletos para o primeiro
                                   % e �ltimo dia na base de dados foram
                                   % retirados previamente.

% Padroniza��o dos dados

% Convers�o da data para o valor do dia do ano correspondente
n = data2dia_vet(Dados(:,1),Dados(:,2),Dados(:,3));

% Corre��o da hora padr�o para o valor da hora solar

B = (n-1)*360/365;
% Equa��o do tempo (em minutos)
E = 229.18*(0.000075 + 0.001868*cosd(B)-0.032077*sind(B)...
           -0.014615*cosd(2*B)-0.04089*sind(2*B));

hora_std = 60*Dados(:,4)+Dados(:,5); % Convers�o da hora padr�o para valor 
                                     % em minutos.
       
hora_solar = hora_std + 4*(long_std - long) + E; % valor total em minutos

omega = -180 + 1/4*hora_solar; % Converte os valores da hora solar em 
                               % minutos para �ngulo hor�rio
                               
d    = 1;       % inicializa��o das variaveis de controle da posi��o dos 
                % elementos nas variaveis de sa�da
hora = 1;

rad_trat(1,1)   = Dados(1,6);   % Valores de radia��o em mJ/m�
n_trat(1)       = n(1);         % Dia do ano
ano_trat(1)     = Dados(1,3);   % Ano do respectivo dia
omega_trat(1,1) = omega(1);     % Valor do �ngulo hor�rio no dia e hora
                                % respectivo
                                
for i = 2:size(Dados,1)   % Para cada elemento nos dados originais
    if n(i) ~= n_trat(d)  % Se o dia da linha � diferente do da linha 
                          % anterior.
        d = d + 1;        % incrementa o �ndice d
        hora = 1;         % retorna o �ndice h para 1 (volta para a hora inicial)
        n_trat(d)   = n(i); % Atualiza o valor do dia 
        ano_trat(d) = Dados(i,3); % Atualiza o valor do ano
    else
        hora = hora + 1; % Se n�o h� mudan�a de dia incrementa a hora
    end
    
    omega_trat(d,hora)  = omega(i);  % omega_trat recebe o pr�ximo valor de �ngulo hor�rio
    rad_trat(d,hora)    = Dados(i,6);% rad_trat recebe o valor de radia��o em mJ/m� 
end

omega_trat = omega_trat(:,20:76); % O limites de coluna para os quais o 
                                  % �ngulo est� entre o nascer e por do sol
                                  % foram obtidos manualmente.
rad_trat = rad_trat(:,20:76);

clear long long_std E hora_std hora_solar omega d hora n Dados i lat
                        % Limpa as variaveis que n�o ser�o utilizadas
                        % posteriormente
                        
save 'dados_Trat'
                        % Salva as variaveis necessarias para a pr�xima
                        % etapa.

%% Etapa 2 - Calculo dos valores de radiacao para uma superficie fixa com 
% inclinacao igual a latitude e orientacao em direcao ao norte

clear;
clc;

load 'dados_Trat'
                        % Carrega os dados salvos na c�lula inicial.
% CONSTANTES

rho_g = 0.2;    % Refletividade do solo
G_sc = 1367;    % Constante solar

% Dados Vi�osa
lat = -20.75;       % [�]
altit = 691.61;     % [m]

n = n_trat;

% �ngulos di�rios da orienta��o do sol
[delta,omega_s] = Geom_Solar_Dia(lat, n_trat);   
                                            % delta = Declina��o solar
                                            % omega_s = �ngulo de por do
                                            % sol


omega_1 = omega_trat(:,1:end-1);% Limites inferiores para os valores de hora
omega_2 = omega_trat(:,2:end);  % Limites superiores para os valores de hora
omega = (omega_1 + omega_2)./2;  % Valor intermedi�rio do �ngulo hor�rio

% C�lculo dos �ngulos de eleva��o solar e de azimute do sol
[theta_z,alpha_s,gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

% RADIA��O SOLAR NO TOPO DA ATMOSFERA

% Radia��o solar di�ria global no topo da atmosfera [MJ/m�]
H_o = Rad_Ext_Dia(lat, n, delta, omega_s);
                            % resultado verificado com os limites
                            % observados no grafico da pg 38 do livro do
                            % Duffie (2013)
                            
B = 2*pi/360*(n-1)*360/365;	% Variavel auxiliar
        
% Radiacao extraterrestre em superficie perpendicular a radiacao,
% conforme Duffie e Beckman (2006) e Iqbal (1983)
G_on = G_sc*(1.000110 + 0.034221*cos(B) + 0.001280*sin(B)...
       + 0.000719*cos(2*B) + 0.000077*sin(2*B));
                            % Resultado verificado com o gr�fico presente
                            % na pg 9 do livro do Duffie (2013)

% Radia��o solar extraterrestre em superficie horizontal
G_o = cosd(theta_z).*(G_on'*ones(1,size(theta_z,2)));
G_o(G_o<0)=0;

% Carrega os valores de radia��o incidente em superficie horizontal medidos
% na estacao meteorologica em Vicosa
G_h = rad_trat(:,2:end);
                            
% �ndice de claridade
k_t = G_h./G_o;
k_t(isnan(k_t)) = 0;    % Retira os valores n�o v�lidos da matriz
k_t(isinf(k_t)) = 0;
k_t(k_t>1)=1;
%k_t=ones(size(k_t,1),size(k_t,2));

% Estimativa das fra��es difusa e direta da radia��o solar hor�ria global
[G_b,G_d] = frac_I(G_h, k_t);

% Dados da superf�cie de refer�ncia (orientada em dire��o ao norte e inclina��o igual a latitude)
beta_sup = abs(lat);
gamma_sup = (lat<0)*180;

% C�LCULO DA RADIA��O EM SUPERF�CIE COM ORIENTA��O ARBITR�RIA

% Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
% sobre uma superf�cie horizontal
R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);

% Radia��o em superf�cie arbitr�ria de acordo com o modelo de Perez et al
G_t_perez = Rad_Inc_Dia_Claro(G_h, G_b, G_d, G_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);
G_t_perez(G_t_perez<0)=0;

%% Algoritmos de rastreamento de tempo fixo, limite de erro beta,
% gamma, ang de incidencia, proposto pelo Aurelio e varia��o da
% proposta do Aur�lio

mods_Rastr = {'fix_Int'};%,'fix_Int_Atr','lim_erro_incid'};
legend_Rastr = {'Algorithm based on d\omega'};%,'Algorithm based on delayed/ d\omega', 'Algorithm based on d\theta'};

d_omega = 3.75;    % Intervalo de �ngulo hor�rio entre os instantes de 
                   % dados disponiveis.
crit_max = 6;      % Valor m�ximo para o crit�rio utilizado
crit_min = 3.75;      % Valor m�nimo para o crit�rio utilizado
crit_inc = 5;      % Valor do incremento para o crit�rio utilizado

% Inicializa as variaveis de ganho e de numero de reorientacoes
g_rastr= zeros(length(mods_Rastr),floor((crit_max-crit_min)/crit_inc)+1);
n_reord_m = zeros(length(mods_Rastr),floor((crit_max-crit_min)/crit_inc)+1);
% C1 = 0.9:-0.2:0.1; % Constantes para a equacao de controle do algoritmo
                     % de controle com base nos angulos de zenith e azimute
% C2 = 0.1:0.2:0.9;

% !!!!!!!!! ATEN��O: O SISTEMA EST� VIFICANDO O POSICIONAMENTO DE 15 EM 15
% MIN, � MAIS INTERESSANTE A IMPLEMENTA��O DO SISTEMA COM VARIA��O DA
% ORIENTA��O EM HOR�RIO MAIS PROPICIO

for j = 1:length(mods_Rastr)
    i = 1;
    for crit = crit_min:crit_inc:crit_max
        [G_t_r,beta_sup,gamma_sup,n_reor_gamma,n_reor_beta] = Rad_Inc_Rast_Real(G_h, G_b, G_d, G_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, mods_Rastr{j}, 'inc_Hora', crit);
        % C�lculo do ganho no valor da radia��o solar anual em
        % rela��o ao sistema com inclina��o igual a latitude
        g_rastr(j,i) = sum(sum(G_t_r*3600*d_omega./15.*(G_t_r(:,:)>0 | isreal(G_t_r(:,:))),2))/sum(sum(G_t_perez*3600*d_omega./15.*(G_t_perez(:,:)>0 | isreal(G_t_perez(:,:))),2));
        % C�lculo do n�mero m�dio de reorienta��es di�rias
        n_reord_m(j,i) = sum(n_reor_gamma + n_reor_beta)/length(n);
        i = i + 1;
    end
end

figure('Color',[1 1 1]);
hold on;
grid on;

plot(sum(G_t_perez,2)./sum(G_h,2))
plot(sum(G_t_r,2)./sum(G_h,2))