function [I_t beta_sup gamma_sup] = Rad_Inc_Rast(I, I_b, I_d, I_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, varargin)

%============ Radia��o Solar Hor�ria Global em sistema solar com rastreamento ======%
% Entradas: I           - Radia��o solar hor�ria global em sup. horizontal          %
%           I_b         - Radia��o solar hor�ria direta                             %
%           I_d         - Radia��o solar hor�ria difusa                             %
%           I_o         - Radia��o solar hor�ria extraterrestre                     %
%           R_b         - Raz�o da radia��o direta em sup incl. em rel a sup. hor.  %
%           rho_g       - Albedo do solo                                            %
%           theta_z     - �ngulo de z�nite [�]                                      %
%           alpha_s     - �ngulo de eleva��o solar [�]                              %
%           gamma_s     - �ngulo de azimute solar [�]                               %
%           omega       - �ngulo hor�rio da metade da hora [�]                      %
%           omega_1     - �ngulo hor�rio do inicio da hora [�]                      %
%           omega_2     - �ngulo hor�rio do t�rmino da hora [�]                     %
%           omega_s     - �ngulo hor�rio de p�r-do-sol [�]                          %
%           lat         - Latitude local [�]                                        %
%           n           - Dia do ano                                                %
%           altit       - Altitude local [m]                                        %
%           delta       - Angulo de declinacao solar [�]                            %
%           varargin{1} - Sele��o do crit�rio de reorienta��o:                      %
%                                   'fix_Int' = intervalos fixos de tempo.          %   
%                                   'lim_Erro_Beta' = Limite fixos para erro em     %
%                                               rela��o a inclina��o ideal.         %
%                                   'lim_Erro_Gamma' = Limite fixos para erro em    %
%                                               rela��o a orientacao ideal.         % 
%           varargin{2} - Determina se a frequ�ncia com que a inclina��o            %
%                         ser� atualizada:                                          %
%                                   'inc_Hor' = Varia��o ao longo do dia.           %
%                                   'inc_Dia' = Varia��o s� no �nicio do dia.       %
%           varargin{3} - Caso ('fix_int') = dT - intervalo de tempo entre          %
%                                                 reorienta��es do sistema.         % 
%                         Caso('lim_Erro_Beta') = dErroBeta - Valor maximo de erro  %
%                                                             em rel a posicao otima%
%                         Caso('lim_Erro_Gamma') = dErroGamma - Valor maximo de erro%
%                                                             em rel a posicao otima%
% Sa�da: I_t - Radia��o solar hor�ria global em superf�cie inclnada [MJ/m�]         %
%        beta_sup - Angulo de inclinacao para cada hora de cada dia dos dados de    %
%                   entrada para o algoritmo de rastreamento selecionado            %    
%        gamma_sup - Angulo de azimute da superficie para cada hora de cada dia dos %
%                   dados de entrada para o algoritmo de rastreamento selecionado   %    
%===================================================================================%

% DETERMINAR UM CASO PADR�O QUANDO N�O FOR ESPECIFICADOS OS PARAMETROS DE
% RASTREAMENTO ???

% registro do n�mero de dias e de horas nos vetores
n_dias = size(alpha_s,1);
n_h = size(alpha_s,2);

% Inicializa��o das matrizes com as orienta��es do sistema para cada dia e
% hora
beta_sup = zeros(n_dias,n_h);
gamma_sup = zeros(n_dias,n_h);

if strcmpi(varargin{1},'fix_Int')   % Sela��o do caso com intervalo fixo
    dT = varargin{3}; % Valor em �ngulo hor�rio do intervalo entre reorienta��es
    if strcmpi(varargin{2},'inc_Dia')   % Sele��o da altera��o da inclina��o da superf�cie apenas uma vez por dia
        for i = 1:n_dias % Para cada dia:
            
            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2)); % Define o valor da inclina��o di�ria como a �tima para o meio dia solar
            
            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega est� entre o nascer e o por-do-sol
            j_old = max(ind(1)-1,1);        % Vari�vel auxiliar para definir o ponto inicial dos c�lculos como o �ndice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como refer�ncia para o c�lculo da passagem do tempo at� a pr�xima reorienta��o
            
            gamma_sup(i,:) = gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orienta��o do sistema at� o omega inicial do rastreamento
            
            for j = ind(1):length(omega) % Para cada hora ap�s a hora inicial do rastreamento
                if (omega(j) - omega_old >= dT)  % Se o intervalo que passou for maior que o definido pelo usu�rio:
                    j_old = j;              % Atualiza o valor de j_old (utilizado para definir qual ser� a orienta��o do sistema)
                    omega_old = omega(j);   % Atualiza o valor de angulo horario de refer�ncia.
                end
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orienta��o do sistema
            end          
        end
    elseif strcmpi(varargin{2},'inc_Hora')  % Sele��o da altera��o da inclina��o da superf�cie ao longo do dia
        for i = 1:n_dias % Para cada dia:
                      
            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega est� entre o nascer e o por-do-sol
            j_old = max(ind(1)-1,1);        % Vari�vel auxiliar para definir o ponto inicial dos c�lculos como o �ndice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como refer�ncia para o c�lculo da passagem do tempo at� a pr�xima reorienta��o
            
            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orienta��o do sistema at� o omega inicial do rastreamento
            
            for j = ind(1):length(omega) % Para cada hora ap�s a hora inicial do rastreamento
                if (omega(j) - omega_old >= dT)  % Se o intervalo que passou for maior que o definido pelo usu�rio:
                    j_old = j;              % Atualiza o valor de j_old (utilizado para definir qual ser� a orienta��o do sistema)
                    omega_old = omega(j);   % Atualiza o valor de angulo horario de refer�ncia.
                end
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orienta��o do sistema
                beta_sup(i,j) = theta_z(i,j_old);
            end          
        end
        
    end
     
elseif strcmpi(varargin{1},'lim_Erro_Beta')  % Rastreamento com base no erro do angulo de inclinacao
    dErroBeta = varargin{3};    % Limite m�ximo do erro do angulo de inclicacao
    
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias % Para cada dia:

            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega est� entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Vari�vel auxiliar para definir o ponto inicial dos c�lculos como o �ndice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como refer�ncia para o c�lculo da passagem do tempo at� a pr�xima reorienta��o

            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = (lat<0)*180;   % Define a orienta��o do sistema at� o omega inicial do rastreamento

            for j = ind(1):length(omega)                
                if ( abs(theta_z(i,j) - theta_z(i,j_old)) >= dErroBeta )  % Se o erro for maior que o definido pelo usu�rio
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual ser� a orienta��o do sistema)
                end

                beta_sup(i,j) = abs(theta_z(i,j_old));
            end
        end
    
    elseif strcmpi(varargin{2},'inc_Hora')
        for i = 1:n_dias % Para cada dia:

            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega est� entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Vari�vel auxiliar para definir o ponto inicial dos c�lculos como o �ndice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como refer�ncia para o c�lculo da passagem do tempo at� a pr�xima reorienta��o

            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orienta��o do sistema at� o omega inicial do rastreamento

            for j = ind(1):length(omega)                
                if ( abs(theta_z(i,j) - theta_z(i,j_old)) >= dErroBeta )  % Se o erro for maior que o definido pelo usu�rio
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual ser� a orienta��o do sistema)
                end

                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orienta��o do sistema                
                beta_sup(i,j) = abs(theta_z(i,j_old));
            end
        end
        
    end    
    
elseif strcmpi(varargin{1},'lim_Erro_Gamma')  % Rastreamento com base no erro do angulo de azimute da superficie
    dErroGamma = varargin{3};    % Limite m�ximo do erro do angulo de inclinacao
    
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias % Para cada dia:

            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2)); % Define o valor da inclina��o di�ria como a �tima para o meio dia solar
            
            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega est� entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Vari�vel auxiliar para definir o ponto inicial dos c�lculos como o �ndice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como refer�ncia para o c�lculo da passagem do tempo at� a pr�xima reorienta��o

            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orienta��o do sistema at� o omega inicial do rastreamento

            for j = ind(1):length(omega)                
                if ( abs(gamma_s(i,j) - gamma_s(i,j_old)) >= dErroGamma )  % Se o erro for maior que o definido pelo usu�rio
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual ser� a orienta��o do sistema)
                end
                
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orienta��o do sistema
                
            end
        end
    
    elseif strcmpi(varargin{2},'inc_Hora')
        for i = 1:n_dias % Para cada dia:

            ind = find(omega > -omega_s(i) & omega < omega_s(i)); % Encontra os indices para os quais o valor de omega est� entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Vari�vel auxiliar para definir o ponto inicial dos c�lculos como o �ndice da hora logo antes do nascer do sol
            omega_old = omega(j_old);       % Define o omega inicial como refer�ncia para o c�lculo da passagem do tempo at� a pr�xima reorienta��o

            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orienta��o do sistema at� o omega inicial do rastreamento

            for j = ind(1):length(omega)                
                if ( abs(gamma_s(i,j) - gamma_s(i,j_old)) >= dErroGamma )  % Se o erro for maior que o definido pelo usu�rio
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual ser� a orienta��o do sistema)
                end

                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orienta��o do sistema                
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

% Raz�o entre a radia��o direta incidente sobre superf�cie inclinada e
% sobre uma superf�cie horizontal
R_b_r = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);
R_b_r(R_b_r<0)=0;

I_t = Rad_Inc(I, I_b, I_d, I_o, R_b_r, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);