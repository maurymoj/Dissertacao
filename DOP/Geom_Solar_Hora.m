function [theta_z,alpha_s,gamma_s] = Geom_Solar_Hora(lat, n, omega, varargin)
%================ �ngulos de eleva��o e de azimute =================%
% Entradas: lat     - latitude local [�]                            %
%           n       - dia do ano                                    %
%           omega   - Hora solar [�]                                %
%           delta   - declina��o solar [�](opcional)                %
% Sa�da: alpha_s - �ngulo eleva��o solar [MJ/m�]                    %
%        gamma_s - �ngulo de azimute solar [MJ/m�]                  %
%===================================================================%

if  size(varargin,2) < 1
    B = (n - 1)*360/365; % Vari�vel auxiliar
    delta = 180/pi* ( 0.006918 - 0.399912*cosd(B) + 0.070257*sind(B)...
                    - 0.006758*cosd(2*B) + 0.000907*sind(2*B) ...
                    - 0.002679*cosd(3*B) + 0.00148*sind(3*B) ); % Declina��o solar
elseif size(varargin,2) < 2
    delta = varargin{1};
end

n_horas = size(omega,2);

if size(omega,1)==1
    theta_z = acosd(  cosd(lat)*cosd(delta)'*cosd(omega)...
                    + sind(lat)*sind(delta)'*ones(1,n_horas) );
            
    alpha_s = 90 - theta_z;

    gamma_s = ( ones(length(n),1)*sign(omega)).*...
        abs(acosd( (cosd(theta_z)*sind(lat) - sind(delta)'*ones(1,n_horas))./ ...
               (sind(theta_z)*cosd(lat)) ));
elseif size(omega,1) > 1
    delta = delta'*ones(1,n_horas);
    
    theta_z = acosd(  cosd(lat)*cosd(delta).*cosd(omega)...
                    + sind(lat)*sind(delta) );
            
    alpha_s = 90 - theta_z;

    gamma_s = (sign(omega)).*...
        abs( acosd( (cosd(theta_z)*sind(lat) - sind(delta))./ ...
               (sind(theta_z)*cosd(lat)) ));
end