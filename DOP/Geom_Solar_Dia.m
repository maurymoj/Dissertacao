function [delta omega_s] = Geom_Solar_Dia(lat, n, varargin)

%================ Ângulos da geometria solar no dia ================%
% Entradas: n        - dia do ano                                   %
%           lat      - latitude local [º]                           %
%           varargin - Modelo para calculo da declinacao, 'Cooper'  %
%                      ou 'Spencer'                                 %
% Saídas: delta   - Declinaçao solar [º]                            %
%         omega_s - Ângulo da hora solar para o por-do-sol, [º]     %
% Obs.: Não indicado para latitudes com módulo superior a 66.5º     %
%===================================================================%

if (~isempty(varargin) && strcmpi(varargin{1}, 'Cooper'))
    % Cálculo da declinação solar utilizando a fórmula de Cooper(1969) -
    % equação menos precisa
    delta = 23.45*sind(360*(284+n)/365);
else
    
    % Variável auxiliar
    B = (n - 1)*360/365;

    % Declinação solar - de Spencer(1971), conforme citado por Iqbal(1983) e
    %                    Duffie e Beckman (2006) - fórmula válida para
    %                    latitudes inferiores a 66.5º e apresenta
    %                    resultados mais precisos em relação à equação de
    %                    Cooper.
    delta = 180/pi* ( 0.006918 - 0.399912*cosd(B) + 0.070257*sind(B)...
                    - 0.006758*cosd(2*B) + 0.000907*sind(2*B) ...
                    - 0.002679*cosd(3*B) + 0.00148*sind(3*B) );
end

% Ângulo horário de por-do-sol - Duffie e Beckman (2006)
omega_s = acosd(-tand(lat)*tand(delta));