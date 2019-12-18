% Script rascunho 3.5

clc;clear;

rastr = 0;  % O sistema possui sistema de rastreamento em dois eixos
    mods_Rastr = {'ideal'};
    %mods_Rastr = {'fix_int','lim_erro_beta','lim_erro_gamma','lim_erro_incid'};%,'prop_Aur','prop_Aur_2'};
    %legend_Rastr = {'Algoritmo com base em d\omega','Algoritmo com base em d\beta','Algoritmo com base em d\gamma','Algoritmo com base em d\theta'};%,'Algoritmo ','prop Aur 2'};
    crit_max = 90;  % Valor máximo para o critério utilizado
    crit_min = 0.2; % Valor mínimo para o critério utilizado
    crit_inc = 3;   % Valor do incremento para o critério utilizado
    C1 = 0.9:-0.2:0.1;
    C2 = 0.1:0.2:0.9;
    legend_Rastr_Aur = {'Algoritmo com base em d\omega','C1 = 0,9 e C2 = 0,1', 'C1 = 0,7 e C2 = 0,3', 'C1 = 0,5 e C2 = 0,5', 'C1 = 0,3 e C2 = 0,7', 'C1 = 0,1 e C2 = 0,9'};
    d_omega = 0.1;    % Incremento nos valores do ângulo da hora solar, [o].

if ~rastr % Parâmetros de entrada para otimização da inclinação e orientação
    beta_max = 90;
    beta_min = 0;
    beta_inc = 5;
    gamma_max = 0;
    gamma_min = 180;
    gamma_inc = -5;
end

cid = 7; % Definição da cidade onde sera feita a simulacao

% VALORES PARA CADA CIDADE, NA ORDEM DE LATITUDE - Somente são válidas
% latitudes com módulo menor ou igual a 66.5º
cidades = {'Manaus - AM','Garanhuns - PE','Brasilia - DF','Belo Horizonte - MG','Campo Grande - MS','Joacaba - SC','Sao Gabriel - RS'};
lat_cid = [-3.10          -8.88             -15.77           -19.80                 -20.43               -27.17         -30.33];
altit_cid = [34.36        869.21            1115.25          937.53                 544.51               525.25         120.6];

% Valor de radiação solar diária incidente em superfície horizontal
% fornecido pelo SWERA em parceria com o INPE, valores em MJ/m²
H_cid = [16.82	16.32	16.68	17.60	16.85	15.08	16.74	19.87	20.23	18.78	20.45	17.14
         21.19	20.03	20.84	20.45	17.07	12.72	13.45	17.48	18.69	22.04	24.37	21.23
         21.94	20.71	19.76	20.04	17.69	14.33	17.37	21.15	23.26	24.44	20.84	18.50
         23.98	21.11	20.66	18.72	15.82	13.59	14.55	18.76	18.60	23.63	17.57	22.57
         22.75	21.30	19.60	20.86	16.52	14.69	14.13	18.19	19.89	23.11	25.95	23.59
         22.44	19.66	19.72	17.32	13.41	11.14	11.06	15.40	16.36	20.78	24.75	25.27
         23.81	21.28	19.53	14.41	12.04	9.23	9.84	14.53   17.54	18.26	24.59	25.60];

% ENTRADAS DO PROGRAMA

% CONSTANTES

rho_g = 0.2;    % Refletividade do solo

G_sc = 1367; % Constante solar

% Dias do ano já convertidos para número total
n = [17 47 75 105 135 162 198 228 258 288 318 344]; % Dias médios de cada mês
n_dias_mes = [31 28 31 30 31 30 31 31 30 31 30 31]; % Número de dias do respectivo mês

% Latitude local [º]
lat = lat_cid(cid);

altit = altit_cid(cid);

% Radiação solar diária global

H = H_cid(cid,:);

% GEOMETRIA DA RADIAÇÃO SOLAR

% Ângulos diários da orientação do sol
[delta omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declinação solar
                                            % omega_s = ângulo de por do
                                            % sol

                                            
%----------------- Cálculos para a situação selecionada ------------------%

if ~rastr
    
    omega_1 = floor(-max(omega_s)/15)*15 : 15 : floor(max(omega_s)/15)*15;
    omega_2 =  ceil(-max(omega_s)/15)*15 : 15 : ceil(max(omega_s)/15)*15;
    omega = (omega_1 + omega_2)/2;

    % Cálculo dos ângulos de elevação solar e de azimute do sol
    [theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

    % RADIAÇÃO SOLAR NO TOPO DA ATMOSFERA

    % Radiação solar diária global no topo da atmosfera
    H_o = Rad_Ext_Dia(lat, n, delta, omega_s);

    % Radiação solar horária global no topo da atmosfera
    I_o = Rad_Ext_Hor(lat, n, omega_1, omega_2, delta, omega_s);

    % CÁLCULO DA RADIAÇÃO SOLAR GLOBAL, DIRETA E DIFUSA EM SUPERFÍCIE
    % HORIZONTAL

    % Cálculo da radiação solar global horária a partir da radiação solar global diária
    [I r_t] = Rad_I_de_H(H, lat, n, omega, delta, omega_s);% AVALIAR PRECISÃO DO CALCULO DE r_t e soma ao longo do dia.

    % Índice de claridade
    k_t = I./I_o;
    k_t(isnan(k_t)) = 0;    % Retira os valores não válidos da matriz
    k_t(isinf(k_t)) = 0;

    % Estimativa das frações difusa e direta da radiação solar horária global
    [I_b I_d] = frac_I(I, k_t);
    
    n_dados = 0;
    m_dados = 0;
      
    for beta_sup = beta_min:beta_inc:beta_max

        n_dados = n_dados + 1;

        for gamma_sup = gamma_min:gamma_inc:gamma_max
    
            m_dados = m_dados +1;

            % CÁLCULO DA RADIAÇÃO EM SUPERFÍCIE COM ORIENTAÇÃO ARBITRÁRIA

            % Razão entre a radiação direta incidente sobre superfície inclinada e
            % sobre uma superfície horizontal
            R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);

            % Radiação em superfície arbitrária de acordo com o modelo de Perez et al
            I_t_perez = Rad_Inc(I, I_b, I_d, I_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);

            H_t_perez = sum(I_t_perez,2);
            
            Z_dados(n_dados, m_dados) = sum(H_t_perez.*n_dias_mes')./3.6;
        end
    
        m_dados = 0;
    end
    
else
        
    clima = 1; % Atmosfera padrão

    omega_1 = floor(-max(omega_s)/15)*15 : d_omega : floor(max(omega_s)/15)*15;
    omega_2 =  ceil(-max(omega_s)/15)*15 : d_omega : ceil(max(omega_s)/15)*15;
    omega = (omega_1 + omega_2)/2;

    % Cálculo dos ângulos de elevação solar e de azimute do sol
    [theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

    % RADIAÇÃO SOLAR NO TOPO DA ATMOSFERA

    % Radiação solar diária global no topo da atmosfera
    H_o = Rad_Ext_Dia(lat, n, delta, omega_s);

    % Radiação solar global no topo da atmosfera

    B = 2*pi/360*(n-1)*360/365;	% Variavel auxiliar

    % CALCULO RADIACAO DE DIA CLARO PARA SUPERFICIE HORIZONTAL

    % Radiacao extraterrestre em superficie perpendicular a radiacao,
    % conforme Duffie e Beckman (2006) e Iqbal (1983)
    G_on = G_sc*(1.000110 + 0.034221*cos(B) + 0.001280*sin(B)...
           + 0.000719*cos(2*B) + 0.000077*sin(2*B));

    % Radiação solar extraterrestre em superficie horizontal
    G_o = cosd(theta_z).*(G_on'*ones(1,size(theta_z,2)));

    % Cálculo da radiação solar global, direta e difusa em superfície
    % horizontal

    [G_h G_b G_d] = Rad_Dia_Claro(theta_z, altit/1000, clima, G_on); % Cálculo da radiação total horizontal e parcelas direta e difusa em dia claro

    % Dados da superfície
    beta_sup = abs(lat);
    gamma_sup = (lat<0)*180;

    % CÁLCULO DA RADIAÇÃO EM SUPERFÍCIE COM ORIENTAÇÃO ARBITRÁRIA

    % Razão entre a radiação direta incidente sobre superfície inclinada e
    % sobre uma superfície horizontal
    R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);

    % Radiação em superfície arbitrária de acordo com o modelo de Perez et al
    G_t_perez = Rad_Inc_Dia_Claro(G_h, G_b, G_d, G_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);

    % Algoritmos de rastreamento de tempo fixo, limite de erro beta,
    % gamma, ang de incidencia, proposto pelo Aurelio e variação da
    % proposta do Aurélio

    ganho_Rastr = zeros(length(mods_Rastr),floor((crit_max-crit_min)/crit_inc)+1);
    n_reord_m = zeros(length(mods_Rastr),floor((crit_max-crit_min)/crit_inc)+1);

    for j = 1:length(mods_Rastr)
        i = 1;
        for crit = crit_min:crit_inc:crit_max
            [G_t_r beta_sup gamma_sup n_reor_gamma n_reor_beta] = Rad_Inc_Rast_Dia_Claro(G_h, G_b, G_d, G_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, mods_Rastr{j}, 'inc_Hora', crit);
            % Cálculo do ganho no valor da radiação solar anual em
            % relação ao sistema com inclinação igual a latitude
            ganho_Rastr(j,i) = sum( sum(G_t_r*3600*d_omega./15,2).*n_dias_mes' )./sum( sum(G_t_perez*3600*d_omega./15,2).*n_dias_mes' );
            % Cálculo do número médio de reorientações diárias
            n_reord_m(j,i) = sum(n_reor_gamma + n_reor_beta)/length(n);
            i = i + 1;
        end
    end

    % Estudo de diferentes algoritmos de rastreamento com diferentes
    % pesos para gamma e beta com base no critério do Aurélio

    ganho_Rastr_Aur_2 = zeros(length(C1),floor((crit_max-crit_min)/crit_inc)+1);
    n_reord_m_Aur_2 = zeros(length(C1),floor((crit_max-crit_min)/crit_inc)+1);

    for j = 1:length(C1)
        i = 1;
        for crit = crit_min:crit_inc:crit_max
            [G_t_r beta_sup gamma_sup n_reor_gamma n_reor_beta] = Rad_Inc_Rast_Dia_Claro(G_h, G_b, G_d, G_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, 'prop_Aur_2', 'inc_Hora', crit, C1(j), C2(j));
            ganho_Rastr_Aur_2(j,i) = sum( sum(G_t_r*3600./(d_omega*15),2).*n_dias_mes' )./sum( sum(G_t_perez*3600./(d_omega*15),2).*n_dias_mes' );
            n_reord_m_Aur_2(j,i) = sum(n_reor_gamma + n_reor_beta)/length(n);
            i = i + 1;
        end
    end
   
end


%---------------------- APRESENTAÇÃO DOS RESULTADOS ----------------------%

if ~rastr
    
    % Abre uma nova figura gerando uma superfície com valores de orientação
    % e inclinação no plano x-y e com o valor de radiação no eixo z (armazenados em Z_dados)
    figure; hold on; grid on;
    surf(gamma_min:gamma_inc:gamma_max,beta_min:beta_inc:beta_max,eta*Z_dados);
    colormap('hot');
    colorbar;
    xlabel('Ângulo de azimute (°)');
    ylabel('Ângulo de inclinação (°)');
    view([90 90]);  % ajuste da orientação da superfície para que se veja o plano x-y (inclinação - orientação)
    axis([gamma_max gamma_min beta_min beta_max]); % Delimita a área do gráfico

elseif rastr    
    figure; 
    hold all;
    grid on;

    for j = 1:length(mods_Rastr)
        plot(n_reord_m(j,:), ganho_Rastr(j,:));
    end

    axis([0 90 1 1.5]);
    legend(legend_Rastr{1:length(mods_Rastr)},'Location', 'SouthEast');
    xlabel('Média anual de reorientações diárias');
    ylabel('Ganho em relação à superfície com \beta = \phi');

    figure;
    hold all;
    grid on;

    plot(n_reord_m(1,:), ganho_Rastr(1,:));

    for j = 1:length(C1)
       plot(n_reord_m_Aur_2(j,:),ganho_Rastr_Aur_2(j,:));
    end

    axis([0 90 1 1.5]);
    legend(legend_Rastr_Aur{1:length(C1)},'Location', 'SouthEast');
    xlabel('Média anual de reorientações diárias');
    ylabel('Ganho em relação à superfície com \beta = \phi');
    
end