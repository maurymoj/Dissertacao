function [G_t,beta_sup,gamma_sup,n_reor_gamma,n_reor_beta] = Rad_Inc_Rast_Real(G, G_b, G_d, G_o, R_b, rho_g, theta_z, alpha_s, gamma_s, omega, omega_1, omega_2, omega_s, lat, n, altit, delta, varargin)
% registro do número de dias e de horas nos vetores
n_dias = size(alpha_s,1);
n_h = size(alpha_s,2);

% Inicialização das matrizes com as orientações do sistema para cada dia e
% hora
beta_sup = zeros(n_dias,n_h);
gamma_sup = zeros(n_dias,n_h);

n_reor_beta = ones(1,n_dias); % seta o número de reorientações diárias do angulo de inclinaçao
n_reor_gamma = ones(1,n_dias); % seta o número de reorientações diárias do angulo de azimute

if strcmpi(varargin{1},'fix_Int')   % Selação do caso com intervalo fixo
    dT = varargin{3}; % Valor em ângulo horário do intervalo entre reorientações
    if strcmpi(varargin{2},'inc_Dia')   % Seleção da alteração da inclinação da superfície apenas uma vez por dia
        for i = 1:n_dias % Para cada dia:
            
            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2)); % Define o valor da inclinação diária como a ótima para o meio dia solar
            
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol
            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação
            
            gamma_sup(i,:) = gamma_sup(i,:) + (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento
            
            for j = ind(1):size(omega,2) % Para cada hora após a hora inicial do rastreamento
                if (omega(i,j) - omega_old >= dT)  % Se o intervalo que passou for maior que o definido pelo usuário:
                    j_old = j;              % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    omega_old = omega(i,j);   % Atualiza o valor de angulo horario de referência.
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
            end          
        end
    elseif strcmpi(varargin{2},'inc_Hora')  % Seleção da alteração da inclinação da superfície ao longo do dia
        for i = 1:n_dias % Para cada dia:
                      
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol
            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação
            
            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento
            
            for j = ind(1):size(omega,2) % Para cada hora após a hora inicial do rastreamento
                if (omega(i,j) - omega_old >= dT)  % Se o intervalo que passou for maior que o definido pelo usuário:
                    j_old = j;              % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    omega_old = omega(i,j);   % Atualiza o valor de angulo horario de referência.
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                    n_reor_beta(i) = n_reor_beta(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
                beta_sup(i,j) = theta_z(i,j_old);
            end          
        end
        
    end
%=========================Em construção===================================%

elseif strcmpi(varargin{1},'fix_Int_Atr')  % Rastreamento com base no intervalo de tempo modificado
    dT = varargin{3}; % Valor em ângulo horário do intervalo entre reorientações
    if strcmpi(varargin{2},'inc_Dia')   % Seleção da alteração da inclinação da superfície apenas uma vez por dia
        for i = 1:n_dias % Para cada dia:
            
            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2)); % Define o valor da inclinação diária como a ótima para o meio dia solar
            
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol
            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação
            
            gamma_sup(i,:) = gamma_sup(i,:) + (omega(i,:) <= omega_old).*gamma_s(i,j_old+1);   % Define a orientação do sistema até o omega inicial do rastreamento
            
            for j = ind(1):size(omega,2) % Para cada hora após a hora inicial do rastreamento
                if (omega(i,j) - omega_old >= dT)  % Se o intervalo que passou for maior que o definido pelo usuário:
                    if (omega(i,j)<0)
                        j_old = j+1;              % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    else
                        j_old = j-1;
                    end
                    omega_old = omega(i,j);   % Atualiza o valor de angulo horario de referência.
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
            end          
        end
    elseif strcmpi(varargin{2},'inc_Hora')  % Seleção da alteração da inclinação da superfície ao longo do dia
        for i = 1:n_dias % Para cada dia:
                      
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol
            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação
            
            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento
            
            for j = ind(1):size(omega,2) % Para cada hora após a hora inicial do rastreamento
                if (omega(i,j) - omega_old >= dT)  % Se o intervalo que passou for maior que o definido pelo usuário:
                    if (omega(i,j)<0)
                        j_old = j+1;              % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    else
                        j_old = j-1;
                    end
                    
                    omega_old = omega(i,j);   % Atualiza o valor de angulo horario de referência.
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                    n_reor_beta(i) = n_reor_beta(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
                beta_sup(i,j) = theta_z(i,j_old);
            end          
        end
        
    end   
%============================================================================%    
elseif strcmpi(varargin{1},'lim_Erro_Beta')  % Rastreamento com base no erro do angulo de inclinacao
    dErroBeta = varargin{3};    % Limite máximo do erro do angulo de inclicacao
    
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias % Para cada dia:

            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = (lat<0)*180;   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):size(omega,2)                
                if ( abs(theta_z(i,j) - theta_z(i,j_old)) >= dErroBeta )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_beta(i) = n_reor_beta(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end

                beta_sup(i,j) = abs(theta_z(i,j_old));
            end
        end
    
    elseif strcmpi(varargin{2},'inc_Hora')
        for i = 1:n_dias % Para cada dia:

            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):size(omega,2)                
                if ( abs(theta_z(i,j) - theta_z(i,j_old)) >= dErroBeta )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                    n_reor_beta(i) = n_reor_beta(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
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
            
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):size(omega,2)                
                if ( abs(gamma_s(i,j) - gamma_s(i,j_old)) >= dErroGamma )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end
                
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
                
            end
        end
    
    elseif strcmpi(varargin{2},'inc_Hora')
        for i = 1:n_dias % Para cada dia:

            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):size(omega,2)                
                if ( abs(gamma_s(i,j) - gamma_s(i,j_old)) >= dErroGamma )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                    n_reor_beta(i) = n_reor_beta(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end

                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema                
                beta_sup(i,j) = abs(theta_z(i,j_old));
            end
        end
        
    end
    
%---------CONSTRUCAO--------------------------------------------------------
elseif strcmpi(varargin{1},'lim_Erro_Incid')
    dErroIncid = varargin{3};    % Limite máximo do erro de acordo com o critério proposto pelo Aurélio
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias % Para cada dia:

            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2)); % Define o valor da inclinação diária como a ótima para o meio dia solar
            
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação
            
            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento
            
            theta = acosd( cosd(theta_z(i,j_old)).*cosd(beta_sup(i,j_old)) ...
                         + sind(theta_z(i,j_old)).*sind(beta_sup(i,j_old)).*cosd(gamma_s(i,j_old) - gamma_sup(i,j_old)) );
            
            for j = ind(1):size(omega,2)                
                     
                if ( theta >= dErroIncid )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end
                
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
                
                theta = acosd( cosd(theta_z(i,j)).*cosd(beta_sup(i,j)) ...
                         + sind(theta_z(i,j)).*sind(beta_sup(i,j)).*cosd(gamma_s(i,j) - gamma_sup(i,j)) );
            end
        end
    elseif strcmpi(varargin{2},'inc_Hora')
        for i = 1:n_dias % Para cada dia:
        
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento
            
            theta = acosd( cosd(theta_z(i,j_old)).*cosd(beta_sup(i,j_old)) ...
                         + sind(theta_z(i,j_old)).*sind(beta_sup(i,j_old)).*cosd(gamma_s(i,j_old) - gamma_sup(i,j_old)) );
            
            for j = ind(1):size(omega,2)
                             
                if ( theta >= dErroIncid )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                    n_reor_beta(i) = n_reor_beta(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end

                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema                
                beta_sup(i,j) = abs(theta_z(i,j_old));
                
                theta = acosd( cosd(theta_z(i,j)).*cosd(beta_sup(i,j)) ...
                         + sind(theta_z(i,j)).*sind(beta_sup(i,j)).*cosd(gamma_s(i,j) - gamma_sup(i,j)) );
            end
        end
         
    end
     
elseif strcmpi(varargin{1},'prop_Aur')
    dErroAur = varargin{3};    % Limite máximo do erro de acordo com o critério proposto pelo Aurélio
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias % Para cada dia:

            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2)); % Define o valor da inclinação diária como a ótima para o meio dia solar
            
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):size(omega,2)                
                if ( sqrt( (gamma_s(i,j) - gamma_s(i,j_old)).^2 + (theta_z(i,j) - theta_z(i,j_old)).^2 ) >= dErroAur )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end
                
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
                
            end
        end
    elseif strcmpi(varargin{2},'inc_Hora')
        for i = 1:n_dias % Para cada dia:
        
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):size(omega,2)                
                if ( sqrt( (gamma_s(i,j) - gamma_s(i,j_old)).^2 + (theta_z(i,j) - theta_z(i,j_old)).^2 ) >= dErroAur )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                    n_reor_beta(i) = n_reor_beta(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end

                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema                
                beta_sup(i,j) = abs(theta_z(i,j_old));
            end
        end
         
    end
elseif strcmpi(varargin{1},'prop_Aur_2')
    dErroAur = varargin{3};    % Limite máximo do erro de acordo com o critério proposto pelo Aurélio
    if size(varargin,2) > 3
        C1 = varargin{4};
        C2 = varargin{5};
    else
        C1 = 0.5;
        C2 = 0.5;
    end
    
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias % Para cada dia:

            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2)); % Define o valor da inclinação diária como a ótima para o meio dia solar
            
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):size(omega,2)                
                if ( sqrt( C1*(gamma_s(i,j) - gamma_s(i,j_old)).^2 + C2*(theta_z(i,j) - theta_z(i,j_old)).^2 ) >= dErroAur )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end
                
                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema
                
            end
        end
    elseif strcmpi(varargin{2},'inc_Hora')
        for i = 1:n_dias % Para cada dia:
        
            [ind_x,ind] = find(omega(i,:) > -omega_s(i) & omega(i,:) < omega_s(i)); % Encontra os indices para os quais o valor de omega está entre o nascer e o por-do-sol

            j_old = max(ind(1)-1,1);        % Variável auxiliar para definir o ponto inicial dos cálculos como o índice da hora logo antes do nascer do sol
            omega_old = omega(i,j_old);       % Define o omega inicial como referência para o cálculo da passagem do tempo até a próxima reorientação

            beta_sup(i,:) = beta_sup(i,:) + (omega(i,:) <= omega_old).*theta_z(i,j_old);
            gamma_sup(i,:) = gamma_sup(i,:) + (omega(i,:) <= omega_old).*gamma_s(i,j_old);   % Define a orientação do sistema até o omega inicial do rastreamento

            for j = ind(1):size(omega,2)                
                if ( sqrt( C1*(gamma_s(i,j) - gamma_s(i,j_old)).^2 + C2*(theta_z(i,j) - theta_z(i,j_old)).^2 ) >= dErroAur )  % Se o erro for maior que o definido pelo usuário
                    j_old = j;                      % Atualiza o valor de j_old (utilizado para definir qual será a orientação do sistema)
                    n_reor_gamma(i) = n_reor_gamma(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                    n_reor_beta(i) = n_reor_beta(i) + 1;  % Atualiza o número de reorientações do ângulo de azimute
                end

                gamma_sup(i,j) = gamma_s(i,j_old);  % Atualiza a orientação do sistema                
                beta_sup(i,j) = abs(theta_z(i,j_old));
            end
        end
         
    end
%---------------------------------------------------------------------------
elseif strcmpi(varargin{1},'ideal') % Rastreamento ideal (sempre perpendicular aos raios solares)
    if strcmpi(varargin{2},'inc_Dia')
        for i = 1:n_dias
            beta_sup(i,:) = theta_z(i,ceil(size(theta_z,2)/2));
        end
        gamma_sup = gamma_s;
        n_reor_gamma = n_h*n_reor_gamma;  % Atualiza o número de reorientações do ângulo de azimute
    elseif strcmpi(varargin{2},'inc_Hora')
        beta_sup = theta_z;
        gamma_sup = gamma_s;
        n_reor_gamma = n_h*n_reor_gamma;  % Atualiza o número de reorientações do ângulo de azimute
        n_reor_beta = n_h*n_reor_beta;  % Atualiza o número de reorientações do ângulo de inclinação
    end
end

% Razão entre a radiação direta incidente sobre superfície inclinada e
% sobre uma superfície horizontal
R_b_r = fator_Geom_Rb(lat, n, beta_sup, gamma_sup, omega_1, omega_2, omega_s, delta);
R_b_r(R_b_r<0)=0;

G_t = Rad_Inc_Dia_Claro(G, G_b, G_d, G_o, R_b_r, beta_sup, rho_g, 'Perez', n, altit, theta_z, gamma_s, gamma_sup);