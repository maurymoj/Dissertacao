% Script rascunho 4.0

clc;clear;

rastr = 0;   % O sistema possui sistema de rastreamento em dois eixos
inc_fix = 1; % Estudo de sistema fixo
inc_var = 0; % Estudo de sistema com inclina��o variada ao longo do ano

if rastr;
    %mods_Rastr = {'ideal'};
    mods_Rastr = {'fix_Int','fix_Int_Atr','lim_erro_incid'};
    legend_Rastr = {'Algorithm based on d\omega','Algorithm based on delayed/ d\omega', 'Algorithm based on d\theta'};
    %mods_Rastr = {'fix_int','lim_erro_beta','lim_erro_gamma','lim_erro_incid'};%,'prop_Aur','prop_Aur_2'};
    %legend_Rastr = {'Algoritmo com base em d\omega','Algoritmo com base em d\beta','Algoritmo com base em d\gamma','Algoritmo com base em d\theta'};%,'Algoritmo ','prop Aur 2'};
    crit_max = 80;  % Valor m�ximo para o crit�rio utilizado
    crit_min = 0.1;   % Valor m�nimo para o crit�rio utilizado
    crit_inc = 3;   % Valor do incremento para o crit�rio utilizado
    C1 = 0.9:-0.2:0.1;
    C2 = 0.1:0.2:0.9;
    legend_Rastr_Aur = {'Algoritmo com base em d\omega','C1 = 0,9 e C2 = 0,1', 'C1 = 0,7 e C2 = 0,3', 'C1 = 0,5 e C2 = 0,5', 'C1 = 0,3 e C2 = 0,7', 'C1 = 0,1 e C2 = 0,9'};
    d_omega = 0.1;    % Incremento nos valores do �ngulo da hora solar, [o].
elseif inc_fix  % Par�metros de entrada para otimiza��o da inclina��o e orienta��o
    beta_max = 90;  % Valor m�ximo de inclina��o
    beta_min = 0;   % Valor m�nimo de inclina��o
    beta_inc = 5;   % Incremento no valor de inclina��o
    gamma_max = 0;  % �dem para azimute
    gamma_min = 180;
    gamma_inc = -5;
elseif inc_var
    beta_min = 0;
    beta_inc_1 = 5;
    beta_inc_2 = 1;
    beta_max = 90;
    n_inc = [1 2 4 12];
end

% % Efici�ncia do m�dulo fotovoltaico
% eta = 0.15;

dia_claro = 1; % Simula��o gerando valores de radia��o em dia claro segundo Hottel (1976)

cid = 3; % Defini��o da cidade onde sera feita a simulacao

% VALORES PARA CADA CIDADE, NA ORDEM DE LATITUDE - Somente s�o v�lidas
% latitudes com m�dulo menor ou igual a 66.5�
cidades = {'Manaus - AM','Garanhuns - PE','Brasilia - DF','Belo Horizonte - MG','Campo Grande - MS','Joacaba - SC','Sao Gabriel - RS'};
lat_cid = [-3.10          -8.88             -15.77           -19.80                 -20.43               -27.17         -30.33];
altit_cid = [34.36        869.21            1115.25          937.53                 544.51               525.25         120.6];

% Valor de radia��o solar di�ria incidente em superf�cie horizontal
% fornecido pelo SWERA em parceria com o INPE, valores em MJ/m�
H_cid = [16.82	16.32	16.68	17.60	16.85	15.08	16.74	19.87	20.23	18.78	20.45	17.14
         21.19	20.03	20.84	20.45	17.07	12.72	13.45	17.48	18.69	22.04	24.37	21.23
         21.94	20.71	19.76	20.04	17.69	14.33	17.37	21.15	23.26	24.44	20.84	18.50
         23.98	21.11	20.66	18.72	15.82	13.59	14.55	18.76	18.60	23.63	17.57	22.57
         22.75	21.30	19.60	20.86	16.52	14.69	14.13	18.19	19.89	23.11	25.95	23.59
         22.44	19.66	19.72	17.32	13.41	11.14	11.06	15.40	16.36	20.78	24.75	25.27
         23.81	21.28	19.53	14.41	12.04	 9.23	 9.84	14.53   17.54	18.26	24.59	25.60];

% ENTRADAS DO PROGRAMA

% CONSTANTES

rho_g = 0.2;    % Refletividade do solo

G_sc = 1367; % Constante solar

% Dias do ano j� convertidos para n�mero total
n = [17 47 75 105 135 162 198 228 258 288 318 344]; % Dias m�dios de cada m�s
n_dias_mes = [31 28 31 30 31 30 31 31 30 31 30 31]; % N�mero de dias do respectivo m�s

% Latitude local [�]
lat = lat_cid(cid);

altit = altit_cid(cid);

% Radia��o solar di�ria global

H = H_cid(cid,:);

% GEOMETRIA DA RADIA��O SOLAR

% �ngulos di�rios da orienta��o do sol
[delta omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declina��o solar
                                            % omega_s = �ngulo de por do
                                            % sol

                                            
%----------------- C�lculos para a situa��o selecionada ------------------%

if inc_fix
    
    omega_1 = floor(-max(omega_s)/15)*15 : 15 : floor(max(omega_s)/15)*15;
    omega_2 =  ceil(-max(omega_s)/15)*15 : 15 :  ceil(max(omega_s)/15)*15;
    omega = (omega_1 + omega_2)/2;

    % C�lculo dos �ngulos de eleva��o solar e de azimute do sol
    [theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

    % RADIA��O SOLAR NO TOPO DA ATMOSFERA

    % Radia��o solar di�ria global no topo da atmosfera
    H_o = Rad_Ext_Dia(lat, n, delta, omega_s);

    % Radia��o solar hor�ria global no topo da atmosfera
    I_o = Rad_Ext_Hor(lat, n, omega_1, omega_2, delta, omega_s);

    % C�LCULO DA RADIA��O SOLAR GLOBAL, DIRETA E DIFUSA EM SUPERF�CIE
    % HORIZONTAL

    % C�lculo da radia��o solar global hor�ria a partir da radia��o solar global di�ria
    [I r_t] = Rad_I_de_H(H, lat, n, omega, delta, omega_s);% AVALIAR PRECIS�O DO CALCULO DE r_t e soma ao longo do dia.

    % �ndice de claridade
    k_t = I./I_o;
    k_t(isnan(k_t)) = 0;    % Retira os valores n�o v�lidos da matriz
    k_t(isinf(k_t)) = 0;

    % Estimativa das fra��es difusa e direta da radia��o solar hor�ria global
    [I_b I_d] = frac_I(I, k_t);
    
    n_dados = 0;
    m_dados = 0;
      
    for beta_sup = beta_min:beta_inc:beta_max

        n_dados = n_dados + 1;

        for gamma_sup = gamma_min:gamma_inc:gamma_max
    
            m_dados = m_dados +1;

            % C�LCULO DA RADIA��O EM SUPERF�CIE COM ORIENTA��O ARBITR�RIA

            % Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
            % sobre uma superf�cie horizontal
            R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);

            % Radia��o em superf�cie arbitr�ria de acordo com o modelo de Perez et al
            I_t_perez = Rad_Inc(I, I_b, I_d, I_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);

            %dados_Graf = [dados_Graf; gamma_sup beta_sup sum(sum(I_t_perez,2))];
            H_t_perez = sum(I_t_perez,2);
            
            Z_dados(n_dados, m_dados) = sum(H_t_perez.*n_dias_mes');
        end
    
        m_dados = 0;
    end

elseif inc_var
    omega_1 = floor(-max(omega_s)/15)*15 : 15 : floor(max(omega_s)/15)*15;
    omega_2 =  ceil(-max(omega_s)/15)*15 : 15 :  ceil(max(omega_s)/15)*15;
    omega = (omega_1 + omega_2)/2;

    % C�lculo dos �ngulos de eleva��o solar e de azimute do sol
    [theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

    % RADIA��O SOLAR NO TOPO DA ATMOSFERA

    % Radia��o solar di�ria global no topo da atmosfera
    H_o = Rad_Ext_Dia(lat, n, delta, omega_s);

    % Radia��o solar hor�ria global no topo da atmosfera
    I_o = Rad_Ext_Hor(lat, n, omega_1, omega_2, delta, omega_s);

    % C�LCULO DA RADIA��O SOLAR GLOBAL, DIRETA E DIFUSA EM SUPERF�CIE
    % HORIZONTAL

    % C�lculo da radia��o solar global hor�ria a partir da radia��o solar global di�ria
    [I r_t] = Rad_I_de_H(H, lat, n, omega, delta, omega_s);% AVALIAR PRECIS�O DO CALCULO DE r_t e soma ao longo do dia.

    % �ndice de claridade
    k_t = I./I_o;
    k_t(isnan(k_t)) = 0;    % Retira os valores n�o v�lidos da matriz
    k_t(isinf(k_t)) = 0;

    % Estimativa das fra��es difusa e direta da radia��o solar hor�ria global
    [I_b I_d] = frac_I(I, k_t);
    
    beta_sup = zeros(1,length(n));
    beta_star = zeros(length(n_inc),length(n));
    gamma_sup = 180*(lat<=0);
    
    I_inc = zeros(length(n), length(omega), length(n_inc)+1);
    H_inc = zeros(length(n),length(n_inc)+1);
    I_inc(:,:,1) = I;
    H_inc(:,1) = sum(I_inc(:,:,1),2);%H;
    
    R_b = zeros(length(n), length(omega));
    I_t_perez = zeros(length(n), length(omega));
    H_t_perez = zeros(length(n),1);
    
    
    for i = 1:length(n_inc)
        if n_inc(i) == 1
            beta_ind = 1:12;
        elseif n_inc(i) == 2
            beta_ind = [1:2 9:12;3:8];
        elseif n_inc(i) == 4
            beta_ind = [1:2 12; 3:5; 6:8; 9:11];
        elseif n_inc(i) == 12
            beta_ind = (1:12)';
        end
        
        for j = 1:size(beta_ind,1)
            rad_total = 0;
            
            for beta = beta_min:beta_inc_1:beta_max
                
                % C�LCULO DA RADIA��O EM SUPERF�CIE COM ORIENTA��O ARBITR�RIA
                
                % Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
                % sobre uma superf�cie horizontal
                
                R_b(beta_ind(j,:),:) = fator_Geom_Rb(lat, n(beta_ind(j,:)), beta, gamma_sup, omega_1, omega_2, omega_s(beta_ind(j,:)), delta(beta_ind(j,:)));

                % Radia��o em superf�cie arbitr�ria de acordo com o modelo de Perez et al
                I_t_perez(beta_ind(j,:),:) = Rad_Inc(I(beta_ind(j,:),:), I_b(beta_ind(j,:),:), I_d(beta_ind(j,:),:), I_o(beta_ind(j,:),:), R_b(beta_ind(j,:),:), beta, rho_g, 'Perez', n(beta_ind(j,:)), altit, theta_z(beta_ind(j,:),:), gamma_s(beta_ind(j,:),:), gamma_sup);
                
                %dados_Graf = [dados_Graf; gamma_sup beta_sup sum(sum(I_t_perez,2))];
                H_t_perez(beta_ind(j,:)) = sum(I_t_perez(beta_ind(j,:),:),2);
               
                total = sum(H_t_perez(beta_ind(j,:)).*n_dias_mes(beta_ind(j,:))');
                
                if total > rad_total
                   rad_total = total;
                   beta_sup(beta_ind(j,:)) = beta;
                   I_inc(beta_ind(j,:),:,i+1) = I_t_perez(beta_ind(j,:),:);
                end
                
            end
            
            for beta = max(beta_sup(beta_ind(j))-10,0):beta_inc_2:min(beta_sup(beta_ind(j))+10,90)

                % C�LCULO DA RADIA��O EM SUPERF�CIE COM ORIENTA��O ARBITR�RIA

                % Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
                % sobre uma superf�cie horizontal
                R_b(beta_ind(j,:),:) = fator_Geom_Rb(lat, n(beta_ind(j,:)), beta, gamma_sup, omega_1, omega_2, omega_s(beta_ind(j,:)), delta(beta_ind(j,:)));

                % Radia��o em superf�cie arbitr�ria de acordo com o modelo de Perez et al
                I_t_perez(beta_ind(j,:),:) = Rad_Inc(I(beta_ind(j,:),:), I_b(beta_ind(j,:),:), I_d(beta_ind(j,:),:), I_o(beta_ind(j,:),:), R_b(beta_ind(j,:),:), beta, rho_g, 'Perez', n(beta_ind(j,:)), altit, theta_z(beta_ind(j,:),:), gamma_s(beta_ind(j,:),:), gamma_sup);

                %dados_Graf = [dados_Graf; gamma_sup beta_sup sum(sum(I_t_perez,2))];
                H_t_perez(beta_ind(j,:)) = sum(I_t_perez(beta_ind(j,:),:),2);
                
                total = sum(H_t_perez(beta_ind(j,:),:).*n_dias_mes(beta_ind(j,:))');
                
                if total > rad_total
                   rad_total = total;
                   beta_sup(beta_ind(j,:)) = beta;
                   I_inc(beta_ind(j,:),:,i+1) = I_t_perez(beta_ind(j,:),:);
                end
                
            end
            
        end
        beta_star(i,:) = beta_sup;
        H_inc(:,i+1) = sum(I_inc(:,:,i+1),2);
    end
    
elseif rastr
    if dia_claro == 1
        
        clima = 1; % Atmosfera padr�o
        
        omega_1 = floor(-max(omega_s)/15)*15 : d_omega : floor(max(omega_s)/15)*15;
        omega_2 =  ceil(-max(omega_s)/15)*15 : d_omega : ceil(max(omega_s)/15)*15;
        omega = (omega_1 + omega_2)/2;

        % C�lculo dos �ngulos de eleva��o solar e de azimute do sol
        [theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

        % RADIA��O SOLAR NO TOPO DA ATMOSFERA

        % Radia��o solar di�ria global no topo da atmosfera
        H_o = Rad_Ext_Dia(lat, n, delta, omega_s);      
    
        % CALCULO RADIACAO DE DIA CLARO PARA SUPERFICIE HORIZONTAL
 
        B = 2*pi/360*(n-1)*360/365;	% Variavel auxiliar
        
        % Radiacao extraterrestre em superficie perpendicular a radiacao,
        % conforme Duffie e Beckman (2006) e Iqbal (1983)
        G_on = G_sc*(1.000110 + 0.034221*cos(B) + 0.001280*sin(B)...
		       + 0.000719*cos(2*B) + 0.000077*sin(2*B));
        
        % Radia��o solar extraterrestre em superficie horizontal
        G_o = cosd(theta_z).*(G_on'*ones(1,size(theta_z,2)));
        
        % C�lculo da radia��o solar global, direta e difusa em superf�cie
        % horizontal para dia claro
        
        [G_h G_b G_d] = Rad_Dia_Claro(theta_z, altit/1000, clima, G_on);
        
        % Dados da superf�cie de refer�ncia (orientada em dire��o ao norte e inclina��o igual a latitude)
        beta_sup = abs(lat);
        gamma_sup = (lat<0)*180;
    
        % C�LCULO DA RADIA��O EM SUPERF�CIE COM ORIENTA��O ARBITR�RIA

        % Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
        % sobre uma superf�cie horizontal
        R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);

        % Radia��o em superf�cie arbitr�ria de acordo com o modelo de Perez et al
        G_t_perez = Rad_Inc_Dia_Claro(G_h, G_b, G_d, G_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);

        % Algoritmos de rastreamento de tempo fixo, limite de erro beta,
        % gamma, ang de incidencia, proposto pelo Aurelio e varia��o da
        % proposta do Aur�lio
        
        ganho_Rastr = zeros(length(mods_Rastr),floor((crit_max-crit_min)/crit_inc)+1);
        n_reord_m = zeros(length(mods_Rastr),floor((crit_max-crit_min)/crit_inc)+1);
        
        for j = 1:length(mods_Rastr)
            i = 1;
            for crit = crit_min:crit_inc:crit_max
                [G_t_r beta_sup gamma_sup n_reor_gamma n_reor_beta] = Rad_Inc_Rast_Dia_Claro(G_h, G_b, G_d, G_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, mods_Rastr{j}, 'inc_Hora', crit);
                % C�lculo do ganho no valor da radia��o solar anual em
                % rela��o ao sistema com inclina��o igual a latitude
                ganho_Rastr(j,i) = sum( sum(G_t_r*3600*d_omega./15,2).*n_dias_mes' )./sum( sum(G_t_perez*3600*d_omega./15,2).*n_dias_mes' );
                % C�lculo do n�mero m�dio de reorienta��es di�rias
                n_reord_m(j,i) = sum(n_reor_gamma + n_reor_beta)/length(n);
                i = i + 1;
            end
        end
        
        % Estudo de diferentes algoritmos de rastreamento com diferentes
        % pesos para gamma e beta com base no crit�rio do Aur�lio
        
%         ganho_Rastr_Aur_2 = zeros(length(C1),floor((crit_max-crit_min)/crit_inc)+1);
%         n_reord_m_Aur_2 = zeros(length(C1),floor((crit_max-crit_min)/crit_inc)+1);
%         
%         for j = 1:length(C1)
%             i = 1;
%             for crit = crit_min:crit_inc:crit_max
%                 [G_t_r beta_sup gamma_sup n_reor_gamma n_reor_beta] = Rad_Inc_Rast_Dia_Claro(G_h, G_b, G_d, G_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, 'prop_Aur_2', 'inc_Hora', crit, C1(j), C2(j));
%                 ganho_Rastr_Aur_2(j,i) = sum( sum(G_t_r*3600./(d_omega*15),2).*n_dias_mes' )./sum( sum(G_t_perez*3600./(d_omega*15),2).*n_dias_mes' );
%                 n_reord_m_Aur_2(j,i) = sum(n_reor_gamma + n_reor_beta)/length(n);
%                 i = i + 1;
%             end
%         end
        
    else
        omega_1 = floor(-max(omega_s)/15)*15 : 15 : floor(max(omega_s)/15)*15;
        omega_2 =  ceil(-max(omega_s)/15)*15 : 15 : ceil(max(omega_s)/15)*15;
        omega = (omega_1 + omega_2)/2;

        % C�lculo dos �ngulos de eleva��o solar e de azimute do sol
        [theta_z alpha_s gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

        % RADIA��O SOLAR NO TOPO DA ATMOSFERA

        % Radia��o solar di�ria global no topo da atmosfera
        H_o = Rad_Ext_Dia(lat, n, delta, omega_s);

        % Radia��o solar hor�ria global no topo da atmosfera
        I_o = Rad_Ext_Hor(lat, n, omega_1, omega_2, delta, omega_s);

        % C�LCULO DA RADIA��O SOLAR GLOBAL, DIRETA E DIFUSA EM SUPERF�CIE
        % HORIZONTAL

        % C�lculo da radia��o solar global hor�ria a partir da radia��o solar global di�ria
        [I r_t] = Rad_I_de_H(H, lat, n, omega, delta, omega_s);% AVALIAR PRECIS�O DO CALCULO DE r_t e soma ao longo do dia.

        % �ndice de claridade
        k_t = I./I_o;
        k_t = k_t./k_t;
        k_t(isnan(k_t)) = 0;    % Retira os valores n�o v�lidos da matriz
        k_t(isinf(k_t)) = 0;

        % Estimativa das fra��es difusa e direta da radia��o solar hor�ria global
        [I_b I_d] = frac_I(I, k_t);
               
        % C�LCULO DA RADIA��O EM SUPERF�CIE COM ORIENTA��O ARBITR�RIA
        
        % Dados da superf�cie
        beta_sup = abs(lat);
        gamma_sup = (lat<0)*180;

        % Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
        % sobre uma superf�cie horizontal
        R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);

        R_b(R_b<0)=0; %INTEGRACAO DE Rb COM HORA CONTENDO POR DO SOL NASCER DO SOL

        % Radia��o em superf�cie arbitr�ria de acordo com o modelo de Perez et al
        I_t_perez = Rad_Inc(I, I_b, I_d, I_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);

        
        [I_t_r beta_sup gamma_sup] = Rad_Inc_Rast(I, I_b, I_d, I_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, mods_Rastr{1}, mods_Rastr{2}, 10);
        ganho_Rastr = sum(sum(I_t_r,2))./sum(sum(I_t_perez,2));
        
    end
    
end


%---------------------- APRESENTA��O DOS RESULTADOS ----------------------%

if inc_fix
    
    % Abre uma nova figura gerando uma superf�cie com valores de orienta��o
    % e inclina��o no plano x-y e com o valor de radia��o no eixo z (armazenados em Z_dados)
    figure; hold on; grid on;
    surf(gamma_min:gamma_inc:gamma_max,beta_min:beta_inc:beta_max,Z_dados);
    colormap('hot');
    colorbar;
    xlabel('�ngulo de azimute (�)','FontName','Times New Roman','FontSize',12);
    ylabel('�ngulo de inclina��o (�)','FontName','Times New Roman','FontSize',12);
    % titulo = strcat('Radia��o total anual - ',cidades{cid},'  (', num2str(abs(lat)), '�S)[kWh/m�]');
    % title(titulo);
    view([90 90]);  % ajuste da orienta��o da superf�cie para que se veja o plano x-y (inclina��o - orienta��o)
    axis([gamma_max gamma_min beta_min beta_max]); % Delimita a �rea do gr�fico

elseif inc_var
    
    
    axis([1 12 floor(min(min(H_inc))) ceil(max(max(H_inc)))]);
	xlabel('')
    
elseif rastr
    if dia_claro == 1
        
        figure; 
        hold all;
        grid on;
        
        for j = 1:length(mods_Rastr)
            plot(n_reord_m(j,:), ganho_Rastr(j,:));
        end
        
        axis([0 90 1 1.5]);
        legend(legend_Rastr{1:length(mods_Rastr)},'Location', 'SouthEast');
        xlabel('M�dia anual de reorienta��es di�rias','FontName','Times New Roman','FontSize',12);
        ylabel('Ganho em rela��o � superf�cie com    \beta = \phi','FontName','Times New Roman','FontSize',12);
        
        figure;
        hold all;
        grid on;
        
  %      plot(n_reord_m(1,:), ganho_Rastr(1,:));
        
%         for j = 1:length(C1)
%            plot(n_reord_m_Aur_2(j,:),ganho_Rastr_Aur_2(j,:));
%         end
%         
%         axis([0 90 1 1.5]);
%         legend(legend_Rastr_Aur{1:length(C1)},'Location', 'SouthEast');
%         xlabel('M�dia anual de reorienta��es di�rias');
%         ylabel('Ganho em rela��o � superf�cie com \beta = \phi');
        
    else
        for i=1:length(n)
            figure; hold on; grid on;
            plot(omega, eta*I_t_perez(i,:),omega, eta*I_t_r(i,:),'r');
            titulo = strcat('Radia��o total anual - ',cidades{cid},'  (', num2str(abs(lat)), '�S), mes ',num2str(i));
            title(titulo);
        end
    end
    
end