function [I_t_conc I_t] = Rad_Inc_Conc(I, I_b, I_d, beta_sup, alpha_s, rho_r1, alpha_1, rho_r2, alpha_2, rho_g, R_b)

%=========== Radia��o Solar Hor�ria Global em superf�cie inclinada =============%
% Entradas: I        - Radia��o solar hor�ria global em sup. horizontal         %
%           I_b      - Radia��o solar hor�ria direta                            %
%           I_d      - Radia��o solar hor�ria difusa                            %
%           beta_sup - �ngulo de inclina��o da superf�cie[�]                    % 
%           alpha_s  - �ngulo de eleva��o solar [�]                             %
%           rho_r1   - Reflet�ncia do refletor                                  %
%           alpha_1  - �ngulo de inclina��o do refletor 1 [�]                   %
%           rho_r2   - Reflet�ncia do refletor 2                                %
%           alpha_2  - �ngulo de inclina��o do refletor 2 [�]                   %
%           rho_g    - Albedo do solo                                           %
% Sa�da: I_t - Radia��o solar hor�ria global em superf�cie inclnada [MJ/m�]     %
%===============================================================================%

% Radia��o direta sobre o coletor
I_dir_col = I_b.*sind(alpha_s + beta_sup);

%I_dir_col
% Radia��o difusa sobre o coletor
I_dif_col = I_d.*(1+cosd(beta_sup))/2;

% Radia��o refletida pelo solo
I_ref_gr = rho_g*I.*(1-cosd(beta_sup))/2;

% Radia��o refletida pelo refletor 1
Qui = beta_sup + 2*alpha_1 - alpha_s;

I_ref_r1 = max(rho_r1*I_b.*sind(Qui),0);
%I_ref_r1 = max(rho_r1*I_b.*sind(Qui).*sind(alpha_s - alpha_1),0);

% Radia��o refletida pelo refletor 
Tau = alpha_s + 2*alpha_2 - beta_sup;

I_ref_r2 = max(rho_r2*I_b.*sind(Tau),0);
%I_ref_r2 = max(rho_r2*I_b.*sind(Tau).*cosd(alpha_s + alpha_2),0);

% Radia��o total incidente sobre o m�dulo
I_t_conc = I_dir_col + I_dif_col + I_ref_gr + I_ref_r1 + I_ref_r2;

I_t = I_dir_col + I_dif_col + I_ref_gr;