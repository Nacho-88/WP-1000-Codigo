
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
%   m_He   -> masa de Helio [kg]
%
% Output:
%   F      -> fuerza gravitatoria [N]
%
% Dependencias (definidas fuera de la función):
%   G, M_T, R_T, M_globo
% ================================================================

function F = F_g_globo(z, m_He)

load constantes.mat G M_T R_T M_globo
 
    F = (G * M_T * (m_He + M_globo)) ./ (R_T + z).^2;
end
