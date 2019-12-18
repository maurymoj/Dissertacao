function [I_t_conc I_t] = Rad_Inc_Conc(I, I_b, I_d, beta_sup, alpha_s, rho_r1, alpha_1, rho_r2, alpha_2, rho_g, R_b)

%=========== Radiação Solar Horária Global em superfície inclinada =============%
% Entradas: I        - Radiação solar horária global em sup. horizontal         %
%           I_b      - Radiação solar horária direta                            %
%           I_d      - Radiação solar horária difusa                            %
%           beta_sup - Ângulo de inclinação da superfície[º]                    % 
%           alpha_s  - Ângulo de elevação solar [º]                             %
%           rho_r1   - Refletência do refletor                                  %
%           alpha_1  - Ângulo de inclinação do refletor 1 [º]                   %
%           rho_r2   - Refletência do refletor 2                                %
%           alpha_2  - Ângulo de inclinação do refletor 2 [º]                   %
%           rho_g    - Albedo do solo                                           %
% Saída: I_t - Radiação solar horária global em superfície inclnada [MJ/m²]     %
%===============================================================================%

% Radiação direta sobre o coletor
I_dir_col = I_b.*sind(alpha_s + beta_sup);

%I_dir_col
% Radiação difusa sobre o coletor
I_dif_col = I_d.*(1+cosd(beta_sup))/2;

% Radiação refletida pelo solo
I_ref_gr = rho_g*I.*(1-cosd(beta_sup))/2;

% Radiação refletida pelo refletor 1
Qui = beta_sup + 2*alpha_1 - alpha_s;

I_ref_r1 = max(rho_r1*I_b.*sind(Qui),0);
%I_ref_r1 = max(rho_r1*I_b.*sind(Qui).*sind(alpha_s - alpha_1),0);

% Radiação refletida pelo refletor 
Tau = alpha_s + 2*alpha_2 - beta_sup;

I_ref_r2 = max(rho_r2*I_b.*sind(Tau),0);
%I_ref_r2 = max(rho_r2*I_b.*sind(Tau).*cosd(alpha_s + alpha_2),0);

% Radiação total incidente sobre o módulo
I_t_conc = I_dir_col + I_dif_col + I_ref_gr + I_ref_r1 + I_ref_r2;

I_t = I_dir_col + I_dif_col + I_ref_gr;