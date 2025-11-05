
% F_g_caja.m
% ---------------------------------------------------------------
% Cálculo de la fuerza gravitatoria del payload (caja)
% 
% Fórmula:
%   F_g(z) = (G * M_T * m_caja) / (R_T + z - r_globo(z) + l_c - L_z/2)^2
%
% Ámbito: válido tanto en ascenso como en descenso
%
% Inputs:
%   z      -> altitud real [m]
%   dz_dt  -> velocidad vertical [m/s] (no se usa, pero se mantiene porque en un futuro se podrá diferenciar entre acenso y descenso)
%
%
% Output:
%   F      -> fuerza gravitatoria [N]
%
% Dependencias (definidas en el main estas ctes y variables):
%   G, M_T, R_T, m_caja, l_c, L_z
%   r_globo(z)
% ================================================================

function F = F_g_caja(z, dz_dt) 
    F = (G * M_T * m_caja) ./ ( (R_T + z - r_globo(z) + l_c - L_z/2).^2 );
end



