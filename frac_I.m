function [I_b I_d] = frac_I(I, k_t, varargin)

%============= Radiação Solar Horária Difusa e direta ==============%
% Entradas: I    - Radiação solar horária global [MJ/m²]            %
%           k_t  - Índice de claridade                              %
%           mod  - Modelo a ser usado par cálculo da rad. difusa  %
%                   (Erbs1982 ou BRL2010) -     (opcional)          %
% Saída: I_b - Radiação solar horária direta [MJ/m²]                %
%        I_d - Radiação solar horária difusa [MJ/m²]                %
%===================================================================%

% % Cálculo da radiação solar horária difusa
% if strcmpi(varargin{1}, 'Erbs1982')
    % Calculo de k_d com base no intervalo de k_t (<= 0.22; 0.22 < < 0.8; >= 0.8)
     k_d = (k_t <= 0.22).*(1 - 0.09*k_t)...
         + (k_t > 0.22 & k_t <= 0.8).*(0.9511 - 0.1604*k_t + 4.388* k_t.^2 ...
                                       - 16.638* k_t.^3 + 12.336*k_t.^4)...
         + (k_t > 0.8).*0.165;
 
     I_d = k_d.*I;
 
% elseif strcmpi(varargin{1}, 'BRL2010')

% falta: HSA, alpha_s, K_t, ksi

%   ksi = (nascer do sol).*k_t+1 +...
%         (nascer < t < por do sol).*(k_t+1 + k_t-1)/2 +...
%         (por do sol).*k_t-1;
%
%   k_d = 1/...
%        (1 + exp(-5.38)+6.63*k_t + 0.006*HSA ...
%           - 0.007*alpha_s + 1.75*K_t + 1.31*ksi)

% else % Default - modelo de Erbs et al. (1982)
%     k_d = (k_t <= 0.22).*(1 - 0.09*k_t)...
%         + (k_t > 0.22 & k_t <= 0.8).*(0.9511 - 0.1604*k_t + 4.388* k_t.^2 ...
%                                       - 16.638* k_t.^3 + 12.336*k_t.^4)...
%         + (k_t > 0.8).*0.168;
% 
%     I_d = k_d.*I;    
% end

% Radiação solar horária direta
I_b = I - I_d;