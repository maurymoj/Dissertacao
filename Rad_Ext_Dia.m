function H_o = Rad_Ext_Dia(lat, n, varargin)

%============ Radiação Solar Diária no topo da atmosfera ===========%
% Entradas: n       - dia do ano                                    %
%           lat     - latitude local [º]                            %
%           delta   - declinação solar [º]          (opcional)      %
%           omega_s - hora solar do por-do-sol [º]  (opcional)      %
%           G_sc    - constante solar [W/m²]        (opcional)      %
% Saída: H_o - Radiação solar diária no topo da atmosfera [MJ/m²]   %
%===================================================================%

if  size(varargin,2) < 1
    G_sc = 1367; % Constante solar [W/m²]
    B = (n - 1)*360/365; % Variável auxiliar
    delta = 180/pi* (006918 - 0.399912*cosd(B) + 0.070257*sind(B)...
                    - 0.006758*cosd(2*B) + 0.000907*sind(2*B) ...
                    - 0.002679*cosd(3*B) + 0.00148*sind(3*B) );% Declinação solar
    omega_s = acosd(-tand(lat)*tand(delta)); % Ângulo horário de por-do-sol
elseif size(varargin,2) < 2
    G_sc = 1367;
    delta = varargin{1};
    omega_s = acosd(-tand(lat)*tand(delta));
elseif size(varargin,2) < 3
    G_sc = 1367;
    delta = varargin{1};
    omega_s = varargin{2};    
elseif size(varargin,2) < 4
    G_sc = varargin{3};
    delta = varargin{1};
    omega_s = varargin{2};
end

H_o =  24*3600/pi*G_sc*(1+ 0.033*cosd(360*n/365)).*...
                       (cosd(lat)*cosd(delta).*sind(omega_s)+...
                        pi*omega_s/180*sind(lat).*sind(delta))/1e6;