function I_o = Rad_Ext_Hor(lat, n, omega_1, omega_2, varargin)

%============ Radia��o Solar Hor�ria no topo da atmosfera ==========%
% Entradas: n       - dia do ano                                    %
%           lat     - latitude local [�]                            %
%           omega_1  - �nicio do intervalo de hora solar [�]        %
%           omega_2  - fim do intervalo de hora solar [�]           %
%           delta   - declina��o solar [�]         (opcional)       %
%           omega_s - hora solar do por-do-sol [�] (opcional)       %
%           G_sc    - constante solar [W/m�]       (opcional)       %
% Sa�da: I_o - Radia��o solar hor�ria no topo da atmosfera [MJ/m�]  %
%===================================================================%

if  size(varargin,2) < 1
    G_sc = 1367;                                                        % Constante solar [W/m�]
    B = (n - 1)*360/365;                                                % Vari�vel auxiliar
    delta = 180/pi* ( 0.006918 - 0.399912*cosd(B) + 0.070257*sind(B)...
                    - 0.006758*cosd(2*B) + 0.000907*sind(2*B) ...
                    - 0.002679*cosd(3*B) + 0.00148*sind(3*B) );         % Declina��o solar
    omega_s = acosd(-tand(lat)*tand(delta));                            % �ngulo hor�rio de por-do-sol
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

n_horas = length(omega_1);

% Se omega_1 e omega_2 est�o dentro do intervalo de -omega_s a omega_s
I_o = ((ones(length(n),1)*omega_2 > -omega_s'*ones(1,n_horas) & ones(length(n),1)*omega_1 > -omega_s'*ones(1,n_horas)) & (ones(length(n),1)*omega_2 < omega_s'*ones(1,n_horas) & ones(length(n),1)*omega_1 < omega_s'*ones(1,n_horas))).*(...
      12*3600/pi*G_sc*...
     (cosd(lat)*((1 + 0.033*cosd(360/365*n)).*cosd(delta))'*(sind(omega_2) - sind(omega_1))...
      +pi/180*sind(lat)*((1 + 0.033*cosd(360/365*n)).*sind(delta))'*(omega_2 - omega_1))./1e6);

% Se omega_2 > -omega_s e omega_1 < -omega_s
I_o = I_o + (ones(length(n),1)*omega_2 > -omega_s'*ones(1,n_horas) & ones(length(n),1)*omega_1 < - omega_s'*ones(1,n_horas)).*(...
      12*3600/pi*G_sc*...
      (cosd(lat)*((1 + 0.033*cosd(360/365*n)).*cosd(delta))'*ones(1,n_horas).*(sind(ones(length(n),1)*omega_2) - sind(-omega_s'*ones(1,n_horas)))...
      +pi/180*sind(lat)*((1 + 0.033*cosd(360/365*n)).*sind(delta))'*ones(1,n_horas).*(ones(length(n),1)*omega_2 - (-omega_s)'*ones(1,n_horas)))./1e6);

% Se omega_1 < omega_s e omega_2 > omega_s
I_o = I_o + (ones(length(n),1)*omega_2 > omega_s'*ones(1,n_horas) & ones(length(n),1)*omega_1 < omega_s'*ones(1,n_horas)).*(...
      12*3600/pi*G_sc*...
      (cosd(lat)*((1 + 0.033*cosd(360/365*n)).*cosd(delta))'*ones(1,n_horas).*(sind(omega_s'*ones(1,n_horas)) - sind(ones(length(n),1)*omega_1))...
      +pi/180*sind(lat)*((1 + 0.033*cosd(360/365*n)).*sind(delta))'*ones(1,n_horas).*(omega_s'*ones(1,n_horas) - ones(length(n),1)*omega_1))./1e6);

I_o(I_o < 0) = 0;