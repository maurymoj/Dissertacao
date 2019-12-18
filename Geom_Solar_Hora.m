function [theta_z,alpha_s,gamma_s] = Geom_Solar_Hora(lat, n, omega, varargin)
%================ Ângulos de elevação e de azimute =================%
% Entradas: lat     - latitude local [º]                            %
%           n       - dia do ano                                    %
%           omega   - Hora solar [º]                                %
%           delta   - declinação solar [º](opcional)                %
% Saída: alpha_s - ângulo elevação solar [MJ/m²]                    %
%        gamma_s - ângulo de azimute solar [MJ/m²]                  %
%===================================================================%

if  size(varargin,2) < 1
    B = (n - 1)*360/365; % Variável auxiliar
    delta = 180/pi* ( 0.006918 - 0.399912*cosd(B) + 0.070257*sind(B)...
                    - 0.006758*cosd(2*B) + 0.000907*sind(2*B) ...
                    - 0.002679*cosd(3*B) + 0.00148*sind(3*B) ); % Declinação solar
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