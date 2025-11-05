
% F_g_globo.m
% ---------------------------------------------------------------
% Cálculo de la fuerza gravitatoria del globo (Helio + látex)
%
% Fórmula:
%   F_g(z) = (G * M_T * (m_He + m_latex)) / (R_T + z)^2
%
% Ámbito: solo durante ascenso
%
% Inputs:
%   z      -> altitud real [m]
%   dz_dt  -> velocidad de ascenso [m/s] 
%   m_He   -> masa de Helio [kg]
%
% Output:
%   F      -> fuerza gravitatoria [N]
%
% Dependencias (definidas fuera de la función):
%   G, M_T, R_T, m_latex
% ================================================================

function F = F_g_globo(z, dz_dt, m_He) 
 
    F = (G * M_T * (m_He + m_latex)) ./ (R_T + z).^2;
end
