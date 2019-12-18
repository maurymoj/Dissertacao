function [G_h,G_b,G_d] = Rad_Dia_Claro(theta_z, A, clima, G_on)

%========== Função para estimativa da radiação em um  dia claro ==========%
% Entradas:                                                               %
%       Theta_z = Ângulo de zênite [°]                                    %
%       A       = Altitude [km]                                           %
%       clima   = Classifica o clima da regiao:
%                  1 - Atmosfera padrão
%                  2 - Clima tropical                                     %
%                  3 - Verao de latitude media                            %
%                  4 - Verao sub artico                                   %
%                  5 - Inverno de latitude media                          %
%       G_on    = Radiação solar extraterrestre em sup. normal a radiação %
%                                                                         %
% Saidas:                                                                 %
%       G_h     = Radiação solar global em sup. horizontal                %
%       G_b     = Radiação solar difusa em sup. horizontal                %
%       G_d     = Radiação solar direta em sup. horizontal                %
%=========================================================================%

% Parametros para calculo da radiacao em dia claro, modelo apresentado por
% Hottel (1976), conforme citado por Duffie e Beckman (2006) e por Rabl (1985).

coeff = [ 1.00 1.00 1.00 1.00                   % Atmosfera padrão
          0.95 0.92 0.98 1.02					% Clima tropical
          0.97 0.96 0.98 1.02					% verao de latitude media
          0.99 0.98 0.99 1.01					% verao sub artico
          1.03 1.04 1.01 1.00 ];                % inverno de latitude media

a_o_star = 0.4237 - 0.00821*(6 - A)^2;

a_1_star = 0.5055 + 0.00595*(6.5 - A)^2;

k_star = 0.2711 + 0.01858*(2.5 - A)^2;

a_o = coeff(clima,1)*a_o_star;					% Calculo dos coeficientes da equacao

a_1 = coeff(clima,3)*a_1_star;

k = coeff(clima,4)*k_star;

altit = A*1000;     % CONFIRMAR

m =(theta_z < 70).* 1./cosd(theta_z)+...
   (theta_z >= 70).* (exp(-0.0001184*altit)./ ...
                     ( cosd(theta_z) + 0.5057*(96.080 - theta_z).^-1.634 ));						% massa de ar
                    % ESTUDAR RESULTADOS NUMEROS COMPLEXOS

tal_b = (theta_z < 90).*(a_o + a_1 * exp( - k ./ cosd(theta_z) ) );	% Transmissividade da radiacao direta

G_b = max((G_on'*ones(1,size(theta_z,2))).*cosd(theta_z).*tal_b,0);				% Radiacao direta em superficie horizontal

tal_d = (theta_z < 90).*(0.2710 - 0.2939*tal_b);					% Transmissividade da radiacao difusa

G_d = max((G_on'*ones(1,size(theta_z,2))).*cosd(theta_z).*tal_d,0);				% Radiacao difusa em sup horizontal

if (G_d == 0 | G_b == 0)
    G_h = 0;
else
    G_h = G_d + G_b;						% Radiacao solar global em sup. horizontal
end