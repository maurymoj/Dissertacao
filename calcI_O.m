%========= Cálculo da radiação solar no topo da atmosfera ==========%
% Entradas: n - dia do ano                                          %
%           lat - latitude local                                    %
%           omega - hora solar                                      %
% Saídas: I_o - Radiação solar horária no topo da atmosfera         % 
%         H_o - Radiação solar diária no topo da atmosfera          %
%===================================================================%
n = 1:365;
lat = 43;
omega_2 = -105:15:120; 
omega_1 = -120:15:105;  % o uso de um só valor (o da média das horas) 
                        % apresenta valor superior ao usando os lims do
                        % intervalo, mas para intervalos de só 1 hra é
                        % satisfatório.
omega = (omega_2 + omega_1)/2;
rho_g = 0.2;

% CONSTANTES
G_sc = 1367; % W/m²

% variável auxiliar
B = (n - 1)*360/365;

% CÁLCULO DOS ÂNGULOS

% Declinação [º]: deslocamento angular do sol em relação 
%                 ao plano do equador. Norte positivo, 
%                 sul negativo [-23.45º, 23.45º]

% forma precisa (+/- 0.035):
delta = 180/pi* ( 0.006918 - 0.399912*cosd(B) + 0.070257*sind(B)...
        - 0.006758*cosd(2*B) + 0.000907*sind(2*B) ...
        - 0.002679*cosd(3*B) + 0.00148*sind(3*B) );

% Ângulo horário do pôr do sol (theta_z = 90º)
omega_s = acosd(-tand(lat)*tand(delta));

% Cálculo do número máximo de divisões das horas após o meio-dia solar
%maxDiv = ceil(max(omega_s)/15);

theta_z = acosd( cosd(lat)*cosd(delta)'*cosd(omega) ...
                 + sind(lat)*sind(delta)'*ones(1,16));%VERIFICAR

% RADIAÇÃO SOLAR HORÁRIA NO TOPO DA ATMOSFERA [MJ/m²]
I_o = ((ones(length(n),1)*omega_2 > -omega_s'*ones(1,16) & ones(length(n),1)*omega_1 > -omega_s'*ones(1,16)) & (ones(length(n),1)*omega_2 < omega_s'*ones(1,16) & ones(length(n),1)*omega_1 < omega_s'*ones(1,16))).*(...
      12*3600/pi*G_sc*...
     (cosd(lat)*((1 + 0.033*cosd(360/365*n)).*cosd(delta))'*(sind(omega_2) - sind(omega_1))...
      +pi/180*sind(lat)*((1 + 0.033*cosd(360/365*n)).*sind(delta))'*(omega_2 - omega_1))/1e6);
I_o = I_o + (ones(length(n),1)*omega_2 > -omega_s'*ones(1,16) & ones(length(n),1)*omega_1 < - omega_s'*ones(1,16)).*(...
      12*3600/pi*G_sc*...
      (cosd(lat)*((1 + 0.033*cosd(360/365*n)).*cosd(delta))'*ones(1,16).*(sind(ones(length(n),1)*omega_2) - sind(omega_s'*ones(1,16)))...
      +pi/180*sind(lat)*((1 + 0.033*cosd(360/365*n)).*sind(delta))'*ones(1,16).*(ones(length(n),1)*omega_2 - omega_s'*ones(1,16)))/1e6);
I_o = I_o + (ones(length(n),1)*omega_2 > omega_s'*ones(1,16) & ones(length(n),1)*omega_1 < omega_s'*ones(1,16)).*(...
      12*3600/pi*G_sc*...
      (cosd(lat)*((1 + 0.033*cosd(360/365*n)).*cosd(delta))'*ones(1,16).*(sind(omega_s'*ones(1,16)) - sind(ones(length(n),1)*omega_1))...
      +pi/180*sind(lat)*((1 + 0.033*cosd(360/365*n)).*sind(delta))'*ones(1,16).*(omega_s'*ones(1,16) - ones(length(n),1)*omega_1))/1e6);
I_o(I_o < 0) = 0;
%for i=1:365
%    I_o(i,omega_2 < -omega_s(i) | omega_1 > omega_s(i),:) = 0;
%end

% RADIAÇÃO SOLAR DIÁRIA NO TOPO DA ATMOSFERA [MJ/m²]
H_o = 24*3600/pi*G_sc*(1+ 0.033*cosd(360*n/365)).*...
                       (cosd(lat)*cosd(delta).*sind(omega_s)+...
                        pi*omega_s/180*sind(lat).*sind(delta))/1e6;

%% RADIAÇÃO SOLAR HORÁRIA GLOBAL EM SUPERFICIE COM INCLINACAO BETA E
%  ORIENTACAO GAMMA

%============= Cálculo da radiação solar horária global ==============%
% Entradas: lat - latitude local   [º]                                %
%           altit - Altitude local [m]                                %
%           n - dia do ano                                            %
%           I_o - Radiação solar horária no topo da atmosfera [MJ/m²] %
%           omega - hora solar                                        %
%           omega_s - hora solar no por-do-sol                        %
%           delta - Ângulo de declinacao                              %
%           I ou H - Radiação solar horária/diária global             %
%           beta - Ângulo de inclinação da superfície                 %
%           gamma - Ângulo de azimute da superfície                   %
% Saídas: I_t - Radiação solar horária em superficie inclinada [MJ/m²]%   
%=====================================================================%
                    
 
% RADIAÇÃO SOLAR HORÁRIA GLOBAL A PARTIR DA RADIAÇÃO SOLAR DIÁRIA GLOBAL

% Cálculo dos coeficientes a e b para a equação de r_t = I/H
%a = 0.409 + 0.5016*sind(omega_s - 60);
%b = 0.6609 - 0.4767*sind(omega_s - 60);

%r_t = pi/24*(a'*ones(1,16)+b'*cosd(omega)).*(ones(length(n),1)*cosd(omega) - cosd(omega_s)'*ones(1,16))./...
%                           (sind(omega_s)'*ones(1,16) - pi*(omega_s/180.*cosd(omega_s))'*ones(1,16));
%r_t(r_t < 0) = 0;%VERIFICAR

%I = r_t.*(H'*ones(1,16));

% RADIAÇÃO SOLAR HORARIA DIFUSA [W/m²]

% Índice de claridade horário k_t
k_t = I./I_o;

% Modelo BRL - Ridley et al. (2010)
%PsiC = (ones(length(n),1)*omega < omega_s'*ones(1,16) & ones(length(n),1)*omega > -omega_s'*ones(1,16))
%       (



% Modelo de Erbs et al. (1982)
k_d = (k_t <= 0.22).*(1 - 0.09*k_t)...
    + (k_t > 0.22 & k_t <= 0.8).*(0.9511 - 0.1604*k_t + 4.388* k_t.^2 ...
                                  - 16.638* k_t.^3 + 12.336*k_t.^4)...
    + (k_t > 0.8).*0.168;

I_d = k_d.*I;


% Radiação solar horária direta
I_b = I - I_d;


% ANGULO DE INCIDENCIA [º]
theta = acosd( sind(lat)*sind(delta)'*ones(1,16)*cosd(beta)... 
               - sind(delta)'*ones(1,16)*cosd(lat)*sind(beta)*cosd(gamma)...
               + cosd(delta)'*cosd(lat)*cosd(beta)*cosd(omega)...
               + cosd(delta)'*sind(lat)*sind(beta)*cosd(gamma)*cosd(omega)...
               + cosd(delta)'*sind(beta)*sind(gamma)*sind(omega) );
           
% RADIAÇÃO INCIDENTE SOBRE SUPERFICIE INCLINADA - MODELO HDKR [MJ/m²]
f = sqrt(I_b./I);
A_i = I_b./I_o;
I_t = max( ( I_b + I_d.*A_i ).*R_b ...
           + I_d.*( 1 - A_i ).*((1 + cosd(beta_1) )/2).* ...
                              (1 + f.*sind(beta_1/2).^3) ...
           + I*rho_g.*(1-cosd(beta_1))/2,0 );
 
% % MODELO DE PEREZ ET AL. (1990) [MJ/m²]
% 
% a = max(cosd(theta),0);
% b = max(cosd(theta_z),cosd(85));
% 
% % Parâmetro de claridade
% eps_c = ((I_d + I_bn)./I_d + 5.535e-6*theta_z.^6)./(1+ 5.535e-6*theta_z.^3);% FALTA CALCULAR I_bn
% 
% % Massa de ar
% m = exp(-0.0001184*altit)./(cosd(theta_z) + 0.5057*(96.080 - theta_z).^-1.634);
% 
% % Parâmetro de brilho - DELTA
% DELTA = m*I_d./I_on;% FALTA CALCULAR I_on
% 
% % Coeficientes de claridade de acordo com o valor de epsilon
% f = [-0.008     0.588   -0.062  -0.060  0.072   -0.022
%       0.130     0.683   -0.151  -0.019  0.066   -0.029
%       0.330     0.487   -0.221   0.055 -0.064   -0.026
%       0.568     0.187   -0.295   0.109 -0.152    0.014
%       0.873    -0.392   -0.362   0.226 -0.462    0.001
%       1.132    -1.237   -0.412   0.288 -0.823    0.056
%       1.060    -1.600   -0.359   0.264 -1.127    0.131
%       0.678    -0.327   -0.250   0.156 -1.377    0.251];
%   
% % Coeficientes de claridade F_1 e F_2
% F_1 = (eps_c >= 1.000 & eps_c < 1.065).*max(0, (f(1,1) + f(1,2)*DELTA + pi*theta_z./180*f(1,3) ))...
%      +(eps_c >= 1.065 & eps_c < 1.230).*max(0, (f(2,1) + f(2,2)*DELTA + pi*theta_z./180*f(2,3) ))...
%      +(eps_c >= 1.230 & eps_c < 1.500).*max(0, (f(3,1) + f(3,2)*DELTA + pi*theta_z./180*f(3,3) ))...
%      +(eps_c >= 1.500 & eps_c < 1.950).*max(0, (f(4,1) + f(4,2)*DELTA + pi*theta_z./180*f(4,3) ))...
%      +(eps_c >= 1.950 & eps_c < 2.800).*max(0, (f(5,1) + f(5,2)*DELTA + pi*theta_z./180*f(5,3) ))...
%      +(eps_c >= 2.800 & eps_c < 4.500).*max(0, (f(6,1) + f(6,2)*DELTA + pi*theta_z./180*f(6,3) ))...
%      +(eps_c >= 4.500 & eps_c < 6.200).*max(0, (f(7,1) + f(7,2)*DELTA + pi*theta_z./180*f(7,3) ))...
%      +(eps_c >= 6.2).*max(0, (f(8,1) + f(8,2)*DELTA + pi*theta_z./180*f(8,3) ));
%  
% F_2 = (eps_c >= 1.000 & eps_c < 1.065).*(f(1,4) + f(1,5)*DELTA + pi*theta_z./180*f(1,6))...
%      +(eps_c >= 1.065 & eps_c < 1.230).*(f(2,4) + f(2,5)*DELTA + pi*theta_z./180*f(2,6))...
%      +(eps_c >= 1.230 & eps_c < 1.500).*(f(3,4) + f(3,5)*DELTA + pi*theta_z./180*f(3,6))...
%      +(eps_c >= 1.500 & eps_c < 1.950).*(f(4,4) + f(4,5)*DELTA + pi*theta_z./180*f(4,6))...
%      +(eps_c >= 1.950 & eps_c < 2.800).*(f(5,4) + f(5,5)*DELTA + pi*theta_z./180*f(5,6))...
%      +(eps_c >= 2.800 & eps_c < 4.500).*(f(6,4) + f(6,5)*DELTA + pi*theta_z./180*f(6,6))...
%      +(eps_c >= 4.500 & eps_c < 6.200).*(f(7,4) + f(7,5)*DELTA + pi*theta_z./180*f(7,6))...
%      +(eps_c >= 6.2).*(f(8,4) + f(8,5)*DELTA + pi*theta_z./180*f(8,6));
% 
% I_t = I_b.*R_b + I_d.*(1 - F_1)*(1+cosd(beta))/2 + I_d.*F_1*a/b...
%       + I_d.*F_2*sind(beta) + I*rho_g*(1-cosd(beta))/2;