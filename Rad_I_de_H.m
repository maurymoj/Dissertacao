function [I r_t] = Rad_I_de_H(H, lat, n, omega, varargin)

%=========== Radia��o Solar Hor�ria Global a partir de H ===========%
% Entradas: H       - Radia��o solar di�ria global [MJ/m�]          %
%           lat     - latitude local [�]                            %
%           n       - dia do ano                                    %
%           omega   - Hora solar [�]                                %
%           delta   - declina��o solar [�](opcional)                %
%           omega_s - hora solar do por-do-sol [�] (opcional)       %
% Sa�da: I - Radia��o solar hor�ria global [MJ/m�]                  %
%===================================================================%

if  size(varargin,2) < 1
    B = (n - 1)*360/365; % Vari�vel auxiliar
    delta = 180/pi* ( 0.006918 - 0.399912*cosd(B) + 0.070257*sind(B)...
                    - 0.006758*cosd(2*B) + 0.000907*sind(2*B) ...
                    - 0.002679*cosd(3*B) + 0.00148*sind(3*B) );% Declina��o solar
    omega_s = acosd(-tand(lat)*tand(delta)); % �ngulo hor�rio de por-do-sol
elseif size(varargin,2) < 2
    delta = varargin{1};
    omega_s = acosd(-tand(lat)*tand(delta));
elseif size(varargin,2) < 3
    omega_s = varargin{2};
end

n_horas = size(omega,2);

% C�lculo dos coeficientes a e b para a equa��o de r_t = I/H
a = 0.409 + 0.5016*sind(omega_s - 60);
b = 0.6609 - 0.4767*sind(omega_s - 60);

r_t = pi/24*(a'*ones(1,n_horas)+b'*cosd(omega)).*(ones(length(n),1)*cosd(omega) - cosd(omega_s)'*ones(1,n_horas))./...
                           (sind(omega_s)'*ones(1,n_horas) - pi*(omega_s/180.*cosd(omega_s))'*ones(1,n_horas));
r_t(r_t < 0) = 0;%VERIFICAR

I = r_t.*(H'*ones(1,n_horas));