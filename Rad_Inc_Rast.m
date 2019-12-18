function [I_t beta_sup gamma_sup] = Rad_Inc_Rast(I, I_b, I_d, I_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, varargin)

%============ Radiação Solar Horária Global em sistema solar com rastreamento ======%
% Entradas: I           - Radiação solar horária global em sup. horizontal          %
%           I_b         - Radiação solar horária direta                             %
%           I_d         - Radiação solar horária difusa                             %
%           I_o         - Radiação solar horária extraterrestre                     %
%           R_b         - Razão da radiação direta em sup incl. em rel a sup. hor.  %
%           rho_g       - Albedo do solo                                            %
%           theta_z     - Ângulo de zênite [°]                                      %
%           alpha_s     - Ângulo de elevação solar [°]                              %
%           gamma_s     - Ângulo de azimute solar [°]                               %
%           omega       - Ângulo horàrio da metade da hora [°]                      %
%           omega_1     - Ângulo horàrio do inicio da hora [°]                      %
%           omega_2     - Ângulo horàrio do término da hora [°]                     %
%           omega_s     - Ângulo horàrio de pôr-do-sol [°]                          %
%           lat         - Latitude local [°]                                        %
%           n           - Dia do ano                                                %
%           altit       - Altitude local [m]                                        %
%           delta       - Angulo de declinacao solar [°]                            %
%           varargin{1} - Seleção do critério de reorientação:                      %
%                                   'fix_Int' = intervalos fixos de tempo.          %   
%                                   'lim_Erro_Beta' = Limite fixos para erro em     %
%                                               relação a inclinação ideal.         %
%                                   'lim_Erro_Gamma' = Limite fixos para erro em    %
%                                               relação a orientacao ideal.         % 
%           varargin{2} - Determina se a frequência com que a inclinação            %
%                         será atualizada:                                          %
%                                   'inc_Hor' = Variação ao longo do dia.           %
%                                   'inc_Dia' = Variação só no ínicio do dia.       %
%           varargin{3} - Caso ('fix_int') = dT - intervalo de tempo entre          %
%                                                 reorientações do sistema.         % 
%                         Caso('lim_Erro_Beta') = dErroBeta - Valor maximo de erro  %
%                                                             em rel a posicao otima%
%                         Caso('lim_Erro_Gamma') = dErroGamma - Valor maximo de erro%
%                                                             em rel a posicao otima%
% Saída: I_t - Radiação solar horária global em superfície inclnada [MJ/m²]         %
%        beta_sup - Angulo de inclinacao para cada hora de cada dia dos dados de    %
%                   entrada para o algoritmo de rastreamento selecionado            %    
%        gamma_sup - Angulo de azimute da superficie para cada hora de cada dia dos %
%                   dados de entrada para o algoritmo de rastreamento selecionado   %    
%===================================================================================%

% DETERMINAR UM CASO PADRÃO QUANDO NÃO FOR ESPECIFICADOS OS PARAMETROS DE
% RASTREAMENTO ???

% registro do número de dias e de horas nos vetores
n_dias = size(alpha_s,1);
n_h = size(alpha_s,2);

% Inicialização das matrizes com as orientações do sistema para cada dia e
% hora
beta_sup = zeros(n_dias,n_h);
gamma_sup = zeros(n_dias,n_h);

if strcmpi(varargin{1},'fix_Int')   % Selação do caso com intervalo fixo
    dT = varargin{3}; % Valor em ângulo horário do intervalo entre reorientações
    if strcmpi(varargin{2},'inc_Dia')   % Seleção da alteração da inclinação da superfície apenas uma vez por dia
        for i = 1:n_dias % Para cada dia:
            
            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2)); % Define o valor da inclinação diária como a ótima para o meio dia solar
            
            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol
            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação
            
            gamma_sup(i,:) = gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento
            
            for j = ind(1):length(omega) % Para cada hora após a hora inicial do rastreamento
                if (omega(j) - omega_old >= dT)  % Se o intervalo que passou for maior que o definido pelo usuário:
                    j_old = j;              % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    omega_old = omega(j);   % Atualiza o valor de angulo horario de referência.
                end
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
            end          
        end
    elseif strcmpi(varargin{2},'inc_Hora')  % Seleção da alteração da inclinação da superfície ao longo do dia
        for i = 1:n_dias % Para cada dia:
                      
            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol
            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação
            
            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento
            
            for j = ind(1):length(omega) % Para cada hora após a hora inicial do rastreamento
                if (omega(j) - omega_old >= dT)  % Se o intervalo que passou for maior que o definido pelo usuário:
                    j_old = j;              % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    omega_old = omega(j);   % Atualiza o valor de angulo horario de referência.
                end
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
                beta_sup(i,j) = theta_z(i,j_old);
            end          
        end
        
    end
     
elseif strcmpi(varargin{1},'lim_Erro_Beta')  % Rastreamento com base no erro do angulo de inclinacao
    dErroBeta = varargin{3};    % Limite máximo do erro do angulo de inclicacao
    
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias % Para cada dia:

            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = (lat<0)*180;   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):length(omega)                
                if ( abs(theta_z(i,j) - theta_z(i,j_old)) >= dErroBeta )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                end

                beta_sup(i,j) = abs(theta_z(i,j_old));
            end
        end
    
    elseif strcmpi(varargin{2},'inc_Hora')
        for i = 1:n_dias % Para cada dia:

            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):length(omega)                
                if ( abs(theta_z(i,j) - theta_z(i,j_old)) >= dErroBeta )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                end

                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema                
                beta_sup(i,j) = abs(theta_z(i,j_old));
            end
        end
        
    end    
    
elseif strcmpi(varargin{1},'lim_Erro_Gamma')  % Rastreamento com base no erro do angulo de azimute da superficie
    dErroGamma = varargin{3};    % Limite máximo do erro do angulo de inclinacao
    
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias % Para cada dia:

            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2)); % Define o valor da inclinação diária como a ótima para o meio dia solar
            
            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):length(omega)                
                if ( abs(gamma_s(i,j) - gamma_s(i,j_old)) >= dErroGamma )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                end
                
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
                
            end
        end
    
    elseif strcmpi(varargin{2},'inc_Hora')
        for i = 1:n_dias % Para cada dia:

            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):length(omega)                
                if ( abs(gamma_s(i,j) - gamma_s(i,j_old)) >= dErroGamma )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                end

                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema                
                beta_sup(i,j) = abs(theta_z(i,j_old));
            end
        end
        
    end
    
%---------CONSTRUCAO--------------------------------------------------------
elseif strcmpi(varargin{1},'prop_??')
%     if strcmpi(varargin{2},'inc_Dia')
%     elseif strcmpi(varargin{2},'inc_Hor')
%     end

%---------------------------------------------------------------------------
elseif strcmpi(varargin{1},'ideal') % Rastreamento ideal (sempre perpendicular aos raios solares)
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias
            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2));
        end
        gamma_sup = gamma_s;
    elseif strcmpi(varargin{2},'inc_Hora')
        beta_sup = theta_z;
        gamma_sup = gamma_s;
    end
end

% Razão entre a radiação direta incidente sobre superfície inclinada e
% sobre uma superfície horizontal
R_b_r = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);
R_b_r(R_b_r<0)=0;

I_t = Rad_Inc(I, I_b, I_d, I_o, R_b_r, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);