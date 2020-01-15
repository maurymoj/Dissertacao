function [delta omega_s] = Geom_Solar_Dia(lat, n, varargin)

%================ �ngulos da geometria solar no dia ================%
% Entradas: n        - dia do ano                                   %
%           lat      - latitude local [�]                           %
%           varargin - Modelo para calculo da declinacao, 'Cooper'  %
%                      ou 'Spencer'                                 %
% Sa�das: delta   - Declina�ao solar [�]                            %
%         omega_s - �ngulo da hora solar para o por-do-sol, [�]     %
% Obs.: N�o indicado para latitudes com m�dulo superior a 66.5�     %
%===================================================================%

if (~isempty(varargin) && strcmpi(varargin{1}, 'Cooper'))
    % C�lculo da declina��o solar utilizando a f�rmula de Cooper(1969) -
    % equa��o menos precisa
    delta = 23.45*sind(360*(284+n)/365);
else
    
    % Vari�vel auxiliar
    B = (n - 1)*360/365;

    % Declina��o solar - de Spencer(1971), conforme citado por Iqbal(1983) e
    %                    Duffie e Beckman (2006) - f�rmula v�lida para
    %                    latitudes inferiores a 66.5� e apresenta
    %                    resultados mais precisos em rela��o � equa��o de
    %                    Cooper.
    delta = 180/pi* ( 0.006918 - 0.399912*cosd(B) + 0.070257*sind(B)...
                    - 0.006758*cosd(2*B) + 0.000907*sind(2*B) ...
                    - 0.002679*cosd(3*B) + 0.00148*sind(3*B) );
end

% �ngulo hor�rio de por-do-sol - Duffie e Beckman (2006)
omega_s = acosd(-tand(lat)*tand(delta));