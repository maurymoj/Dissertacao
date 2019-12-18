function I_t = Rad_Inc(I, I_b, I_d, I_o, R_b, beta_sup, rho_g, mod, varargin)

%============ Radiação Solar Horária Global em superfície inclinada ===============%
% Entradas: I           - Radiação solar horária global em sup. horizontal         %
%           I_b         - Radiação solar horária direta                            %
%           I_d         - Radiação solar horária difusa                            %
%           I_o         - Radiação solar horária extraterrestre                    %
%           R_b         - Razão da radiação direta em sup incl. em rel a sup. hor. %
%           beta_sup    - Ângulo de inclinação da superfície [º]                   %
%           rho_g       - Albedo do solo                                           %
%           mod         - Modelo de radiação em sup. inc.('HDKR' ou 'Perez')       %
%           varargin{1} - para Perez -> dia do ano                                 %
%           varargin{2} - para Perez -> Altitude local                             %
%           varargin{3} - para Perez -> Ângulo de zênite                           %
%           +                                                                      %
%           varargin{4} - para Peraz -> Ângulo de incidência                       %
%           ou                                                                     %
%           varargin{4} - para Perez -> Ângulo de azimute solar                    %
%           varargin{5} - para Perez -> Ângulo de azimute da superficie            %
%                                                                                  %
% Saída: I_t - Radiação solar horária global em superfície inclnada [MJ/m²]        %
%==================================================================================%

G_sc = 1367; % Constante solar  - UNIVERSALISAR USO DA CONSTANTE SOLAR

if strcmpi(mod, 'HDKR')
    % Modelo HDKR - não aconselhavel para gamma (ângulo de azimute da 
    % superficie) distante de 0º para o hem. norte ou 180º para o 
    % hemisfério sul.
    f = sqrt(I_b./I);
    A_i = I_b./I_o;
    I_t = max( ( I_b + I_d.*A_i ).*R_b ...
               + I_d.*( 1 - A_i ).*((1 + cosd(beta_sup) )/2).* ...
                                  (1 + f.*sind(beta_sup/2).^3) ...
               + I*rho_g.*(1-cosd(beta_sup))/2,0 );
elseif strcmpi(mod, 'Perez')
    % Modelo de Perez et al (1990) - Mais preciso mas mais complexo e
    % necessita de um número maior de variaveis de entrada.
    
    n_horas = size(I,2);
    
    if size(varargin,2) < 5
        n = varargin{1};
        altit = varargin{2};
        theta_z = varargin{3};
        theta = varargin{4};
        %n_dias = length(n);
    elseif size(varargin,2) < 6
        n = varargin{1};
        altit = varargin{2};
        theta_z = varargin{3};
        gamma_s = varargin{4};
        gamma_sup = varargin{5};
        theta = acosd( cosd(theta_z).*cosd(beta_sup) ...
                     + sind(theta_z).*sind(beta_sup).*cosd(gamma_s - gamma_sup) );
    end
    
    theta = real(theta);        % VERIFICAR
    %altit = varargin{2};
    a = max(cosd(theta),0);
    b = max(cosd(theta_z),cosd(85));

    I_bn = I_b./cosd(theta_z);         %CONFIRMAR USO
    
    % Parâmetro de claridade
    eps_c = ((I_d + I_bn)./I_d + 5.535e-6*theta_z.^6)./(1+ 5.535e-6*theta_z.^3);

    % Massa de ar
    m = (theta_z < 70).* 1./cosd(theta_z)+...
        (theta_z >= 70).* (exp(-0.0001184*altit)./ ...
                          ( cosd(theta_z) + 0.5057*(96.080 - theta_z).^-1.634 )); % Estudar resultados em num complexo
    
    % Cálculo da rad extraterrestre incidente em plano normal aos raios
    % solares
    
    B = (n-1)'*360/365*ones(1,n_horas);
    
    I_on = 3600*G_sc/1e6*(1.000110 + 0.034221*cosd(B) + 0.001280*sind(B)+...
                      0.000719*cosd(2*B) + 0.000077*sind(2*B));
    
    % Parâmetro de brilho - DELTA
    DELTA = m.*I_d./I_on;

    % Coeficientes de claridade de acordo com o valor de epsilon
    f = [-0.008     0.588   -0.062  -0.060  0.072   -0.022
          0.130     0.683   -0.151  -0.019  0.066   -0.029
          0.330     0.487   -0.221   0.055 -0.064   -0.026
          0.568     0.187   -0.295   0.109 -0.152    0.014
          0.873    -0.392   -0.362   0.226 -0.462    0.001
          1.132    -1.237   -0.412   0.288 -0.823    0.056
          1.060    -1.600   -0.359   0.264 -1.127    0.131
          0.678    -0.327   -0.250   0.156 -1.377    0.251];

    % Coeficientes de claridade F_1 e F_2
    F_1 = (eps_c >= 1.000 & eps_c < 1.065).*max(0, (f(1,1) + f(1,2)*DELTA + pi*theta_z./180*f(1,3) ))...
         +(eps_c >= 1.065 & eps_c < 1.230).*max(0, (f(2,1) + f(2,2)*DELTA + pi*theta_z./180*f(2,3) ))...
         +(eps_c >= 1.230 & eps_c < 1.500).*max(0, (f(3,1) + f(3,2)*DELTA + pi*theta_z./180*f(3,3) ))...
         +(eps_c >= 1.500 & eps_c < 1.950).*max(0, (f(4,1) + f(4,2)*DELTA + pi*theta_z./180*f(4,3) ))...
         +(eps_c >= 1.950 & eps_c < 2.800).*max(0, (f(5,1) + f(5,2)*DELTA + pi*theta_z./180*f(5,3) ))...
         +(eps_c >= 2.800 & eps_c < 4.500).*max(0, (f(6,1) + f(6,2)*DELTA + pi*theta_z./180*f(6,3) ))...
         +(eps_c >= 4.500 & eps_c < 6.200).*max(0, (f(7,1) + f(7,2)*DELTA + pi*theta_z./180*f(7,3) ))...
         +(eps_c >= 6.2).*max(0, (f(8,1) + f(8,2)*DELTA + pi*theta_z./180*f(8,3) ));

    F_2 = (eps_c >= 1.000 & eps_c < 1.065).*(f(1,4) + f(1,5)*DELTA + pi*theta_z./180*f(1,6))...
         +(eps_c >= 1.065 & eps_c < 1.230).*(f(2,4) + f(2,5)*DELTA + pi*theta_z./180*f(2,6))...
         +(eps_c >= 1.230 & eps_c < 1.500).*(f(3,4) + f(3,5)*DELTA + pi*theta_z./180*f(3,6))...
         +(eps_c >= 1.500 & eps_c < 1.950).*(f(4,4) + f(4,5)*DELTA + pi*theta_z./180*f(4,6))...
         +(eps_c >= 1.950 & eps_c < 2.800).*(f(5,4) + f(5,5)*DELTA + pi*theta_z./180*f(5,6))...
         +(eps_c >= 2.800 & eps_c < 4.500).*(f(6,4) + f(6,5)*DELTA + pi*theta_z./180*f(6,6))...
         +(eps_c >= 4.500 & eps_c < 6.200).*(f(7,4) + f(7,5)*DELTA + pi*theta_z./180*f(7,6))...
         +(eps_c >= 6.2).*(f(8,4) + f(8,5)*DELTA + pi*theta_z./180*f(8,6));

    % Cálculo da radiação spobre a superfície
    I_t = I_b.*R_b + I_d.*(1 - F_1).*(1+cosd(beta_sup))/2 + I_d.*F_1.*a./b...
          + I_d.*F_2.*sind(beta_sup) + I*rho_g.*(1-cosd(beta_sup))/2;
else % Modelo padrão - modelo HDKR
    f = sqrt(I_b./I);
    A_i = I_b./I_o;
    I_t = max( ( I_b + I_d.*A_i ).*R_b ...
               + I_d.*( 1 - A_i ).*((1 + cosd(beta_sup) )/2).* ...
                                  (1 + f.*sind(beta_sup/2).^3) ...
               + I*rho_g.*(1-cosd(beta_sup))/2,0 );
end