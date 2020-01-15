clear;
clc;

mods_Rastr = {'ideal','fix_Int_BiAxial_DO_1','fix_Int_BiAxial_DO_2','fix_Int_UniAxial_DO_1','fix_Int_UniAxial_DO_2'};

n_reor = [3:9]; % N�mero de reorienta��es
% n_reor = [3];

crit_max = 80;  % Valor m�ximo para o crit�rio utilizado
crit_min = 0.1;   % Valor m�nimo para o crit�rio utilizado
crit_inc = 3;   % Valor do incremento para o crit�rio utilizado

d_omega = 1;    % Incremento nos valores do �ngulo da hora solar, [o].
dia_claro = 1; % Se 1, Simula��o gerando valores de radia��o em dia claro segundo Hottel (1976)
	           % Se 0, dados reais

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
rho_g = 0.2; % Refletividade do solo
G_sc = 1367; % Constante solar

% Dias do ano j� convertidos para n�mero total
n = [17 47 75 105 135 162 198 228 258 288 318 344]; % Dias m�dios de cada m�s
n_dias_mes = [31 28 31 30 31 30 31 31 30 31 30 31]; % N�mero de dias do respectivo m�s


for cid = 1%:3:7
    for n_r=1:length(n_reor)
        crit = [1,n_reor(n_r),n_reor(n_r),n_reor(n_r),n_reor(n_r)];

    %     hAx = figure('color',[1 1 1]);
%         hA = 1;
    
    %                cid = 7; % Defini��o da cidade onde sera feita a simulacao
        % Latitude local [�]
        lat = lat_cid(cid);
        altit = altit_cid(cid);

        % Radia��o solar di�ria global
        H = H_cid(cid,:);

        % GEOMETRIA DA RADIA��O SOLAR

        % �ngulos di�rios da orienta��o do sol
        [delta,omega_s] = Geom_Solar_Dia(lat, n);   % delta = Declina��o solar
                                                    % omega_s = �ngulo de por do
                                                    % sol

        if dia_claro == 1
            clima = 1; % Atmosfera padr�o

            omega_1 = floor(-max(omega_s)/15)*15 : d_omega : floor(max(omega_s)/15)*15;
            omega_2 =  ceil(-max(omega_s)/15)*15 : d_omega : ceil(max(omega_s)/15)*15;
            omega = (omega_1 + omega_2)/2;

            % C�lculo dos �ngulos de eleva��o solar e de azimute do sol
            [theta_z,alpha_s,gamma_s] = Geom_Solar_Hora(lat, n, omega, delta);

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

            [G_h,G_b,G_d] = Rad_Dia_Claro(theta_z, altit/1000, clima, G_on);

            % Dados da superf�cie de refer�ncia (orientada em dire��o ao norte e inclina��o igual a latitude)
            beta_sup = abs(lat);
            
            gamma_sup = sign(gamma_s).*((lat<0).*180*ones(size(gamma_s,1),size(gamma_s,2)));

            % C�LCULO DA RADIA��O EM SUPERF�CIE COM ORIENTA��O ARBITR�RIA

            % Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
            % sobre uma superf�cie horizontal
            R_b = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);

            % Radia��o em superf�cie arbitr�ria de acordo com o modelo de Perez et al
            G_t_perez = Rad_Inc_Dia_Claro(G_h, G_b, G_d, G_o, R_b, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);
            G_t_rastr = zeros(size(G_t_perez,1),size(G_t_perez,2),length(mods_Rastr));

            
%======================= Rastreamento ==============================                       
            for j=1:length(mods_Rastr)

                [G_t_rastr(:,:,j),beta_track(:,:,j),gamma_track(:,:,j),n_reor_gamma,n_reor_beta] = Rad_Inc_Rast_Dia_Claro(G_h, G_b, G_d, G_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, mods_Rastr{j}, 'inc_Hora', crit(j));

                % C�lculo do ganho no valor da radia��o solar anual em
                % rela��o ao sistema com inclina��o igual a latitude
                ganho_Rastr(n_r,j) = sum( sum(G_t_rastr(:,:,j)*3600*d_omega./15,2).*n_dias_mes' )./sum( sum(G_t_perez*3600*d_omega./15,2).*n_dias_mes' );
                % C�lculo do n�mero m�dio de reorienta��es di�rias
                n_reord_m(n_r,j) = sum(n_reor_gamma + n_reor_beta)/length(n);
            end

%===================================================================
            
        end

        figure('color',[1 1 1])
        mes = {'March','June','September','December'};
        k=1;
        for i=3:3:12
            subplot(2,2,k)
            plot(omega,G_t_perez(i,:))
            hold on
            grid on
            for j=1:length(mods_Rastr)
                plot(omega,G_t_rastr(i,:,j)) 
                g_Rastr(k,j) = (sum(G_t_rastr(i,:,j)*3600*d_omega./15,2).*n_dias_mes(i))./(sum(G_t_perez(i,:)*3600*d_omega./15,2).*n_dias_mes(i));
            end

            title(mes{k});
            k = k+1;
        end

    end
    
    figure('color',[1 1 1])
    grid on
    hold on
    for j=1:length(mods_Rastr)
        plot(n_reor,ganho_Rastr(:,j))
    end
    legend('Ideal','Relogio Biax','Proposto Biax','Relogio Uni','Proposto Uni')
end